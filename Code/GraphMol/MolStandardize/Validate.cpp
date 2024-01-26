//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Validate.h"
#include "Fragment.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogParams.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/PeriodicTable.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

std::vector<ValidationErrorInfo> CompositeValidation::validate(
  const ROMol &mol, bool reportAllFailures) const
{
  std::vector<ValidationErrorInfo> errors;
  for (const auto & method : validations) {
    auto partial = method->validate(mol, reportAllFailures);
    if (!partial.empty()) {
      std::copy(partial.begin(), partial.end(), std::back_inserter(errors));
      if (!reportAllFailures) {
        break;
      }
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> RDKitValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  ROMol molCopy = mol;
  std::vector<ValidationErrorInfo> errors;

  unsigned int na = mol.getNumAtoms();

  if (!na) {
    errors.emplace_back("ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  // loop over atoms
  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    Atom *atom = molCopy.getAtomWithIdx(i);
    try {
      atom->calcExplicitValence();
    } catch (const MolSanitizeException &e) {
      errors.push_back("INFO: [ValenceValidation] " + std::string(e.what()));
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo>
NoAtomValidation::validate(const ROMol &mol, bool /*reportAllFailures*/) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();
  if (!na) {
    errors.emplace_back("ERROR: [NoAtomValidation] Molecule has no atoms");
  }
  return errors;
}

std::vector<ValidationErrorInfo>
FragmentValidation::validate(const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  // REVIEW: reportAllFailures is not being used here. is that correct?
  RDUNUSED_PARAM(reportAllFailures);
  std::shared_ptr<FragmentCatalogParams> fparams(new FragmentCatalogParams(""));
  FragmentCatalog fcat(fparams.get());

  const std::vector<std::shared_ptr<ROMol>> &fgrps = fparams->getFuncGroups();
  INT_VECT mapping;
  VECT_INT_VECT atom_mapping;
  std::vector<ROMOL_SPTR> frags =
      MolOps::getMolFrags(mol, true, &mapping, &atom_mapping);

  for (auto &fgrp : fgrps) {
    std::string fname;
    fgrp->getProp(common_properties::_Name, fname);
    std::vector<RDKit::MatchVectType> res;
    unsigned int matches = SubstructMatch(mol, *fgrp, res);
    //		std::cout << fname << " matches " << matches << std::endl;
    if (matches != 0 && frags.size() != 0) {
      VECT_INT_VECT substructmap;  // store idxs of frag from substructmatch
      for (const auto &match : res) {
        std::vector<int> vec;
        for (const auto &pair : match) {
          vec.push_back(pair.second);
        }
        substructmap.push_back(vec);
      }

      // to stop the same fragment being reported many times if present
      // multiple times in molecule
      bool fpresent = false;

      for (auto &molfragidx : atom_mapping) {
        std::sort(molfragidx.begin(), molfragidx.end());
        for (auto &substructidx : substructmap) {
          std::sort(substructidx.begin(), substructidx.end());
          //					// help to debug...
          //					std::cout << "molfragidx: "  <<
          // std::endl; 					for (const auto
          // &i : molfragidx)
          // {
          // std::cout << i; }; std::cout
          // << std::endl; std::cout <<
          //"substructidx: "  << std::endl;
          // for (const auto &i : substructidx) { std::cout << i; };
          // std::cout <<
          // std::endl;
          //					//
          if ((molfragidx == substructidx) && !fpresent) {
            std::string msg = fname + " is present";
            errors.push_back("INFO: [FragmentValidation] " + msg);
            fpresent = true;
          }
        }
      }
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo>
NeutralValidation::validate(const ROMol &mol, bool /*reportAllFailures*/) const {
  std::vector<ValidationErrorInfo> errors;
  int charge = RDKit::MolOps::getFormalCharge(mol);
  if (charge != 0) {
    std::string charge_str;
    if (charge > 0) {
      charge_str = "+" + std::to_string(charge);
    } else {
      charge_str = std::to_string(charge);
    }
    std::string msg = "Not an overall neutral system (" + charge_str + ')';
    errors.push_back("INFO: [NeutralValidation] " + msg);
  }
  return errors;
}

std::vector<ValidationErrorInfo>
IsotopeValidation::validate(const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  for (auto atom : mol.atoms()) {
    unsigned int isotope = atom->getIsotope();
    if (isotope == 0) {
      continue;
    }

    std::string symbol = atom->getSymbol();
    if (strict) {
      PeriodicTable *periodicTable = PeriodicTable::getTable();
      double mass = periodicTable->getMassForIsotope(symbol, isotope);
      if (mass == 0.0) {
        errors.push_back(
          "ERROR: [IsotopeValidation] Molecule contains unknown isotope " +
          std::to_string(isotope) + symbol);
      }
    } else {
      errors.push_back("INFO: [IsotopeValidation] Molecule contains isotope " +
                       std::to_string(isotope) + symbol);
    }

    if (!errors.empty() && !reportAllFailures) {
      break;
    }
  }
  return errors;
}

// constructor
MolVSValidation::MolVSValidation()
  : CompositeValidation({
      std::make_shared<NoAtomValidation>(),
      std::make_shared<FragmentValidation>(),
      std::make_shared<NeutralValidation>(),
      std::make_shared<IsotopeValidation>()
      })
{
}

// overloaded constructor
MolVSValidation::MolVSValidation(
    const std::vector<std::shared_ptr<ValidationMethod>> & validations)
  : CompositeValidation(validations)
{
}

std::vector<ValidationErrorInfo> AllowedAtomsValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();

  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *qatom = mol.getAtomWithIdx(i);
    bool match = false;
    // checks to see qatom matches one of list of allowedAtoms
    for (const auto &allowedAtom : this->d_allowedList) {
      if (allowedAtom->Match(qatom)) {
        match = true;
        break;
      }
    }
    // if no match, append to list of errors.
    if (!match) {
      std::string symbol = qatom->getSymbol();
      errors.push_back("INFO: [AllowedAtomsValidation] Atom " + symbol +
                       " is not in allowedAtoms list");
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> DisallowedAtomsValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();

  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *qatom = mol.getAtomWithIdx(i);
    bool match = false;
    // checks to see qatom matches one of list of allowedAtoms
    for (const auto &disallowedAtom : this->d_disallowedList) {
      if (disallowedAtom->Match(qatom)) {
        match = true;
      }
    }
    // if no match, append to list of errors.
    if (match) {
      std::string symbol = qatom->getSymbol();
      errors.push_back("INFO: [DisallowedAtomsValidation] Atom " + symbol +
                       " is in disallowedAtoms list");
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> FeaturesValidation::validate(
      const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  // disallow query and dummy atoms, and aliases
  for (auto atom : mol.atoms()) {
    if (atom->hasQuery()) {
      errors.push_back(
        "ERROR: [FeaturesValidation] Query atom " + std::to_string(atom->getIdx()+1) + " is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }
    else if (atom->getAtomicNum() == 0) { // (atom->isDummyAtom()) { // GH PR # 6768
      errors.push_back(
        "ERROR: [FeaturesValidation] Dummy atom " + std::to_string(atom->getIdx()+1) + " is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }

    if (atom->hasProp(common_properties::molFileAlias)) {
      errors.push_back(
        "ERROR: [FeaturesValidation] Atom " + std::to_string(atom->getIdx()+1)
        + " with alias '" + atom->getProp<std::string>(common_properties::molFileAlias)
        + "' is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }
  }

  // disallow query bonds
  for (auto bond : mol.bonds()) {
    if (bond->hasQuery()) {
      errors.push_back(
        "ERROR: [FeaturesValidation] Query bond " + std::to_string(bond->getIdx()+1) + " is not allowed");
      if (!reportAllFailures) {
        return errors;
      }
    }
  }

  // disallow using the enahanced stereochemistry
  // (based on configured options)
  if (!allowEnhancedStereo && mol.getStereoGroups().size()) {
    errors.emplace_back(
      "ERROR: [FeaturesValidation] Enhanced stereochemistry features are not allowed");
  }

  return errors;
}

std::vector<ValidationErrorInfo> Is2DValidation::validate(
      const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  if (!mol.getNumConformers()) {
    errors.emplace_back("ERROR: [Is2DValidation] Molecule has no conformer");
    return errors;
  }

  const auto & conf = mol.getConformer();

  if (conf.is3D()) {
    errors.emplace_back("ERROR: [Is2DValidation] Molecule is 3D");
    return errors;
  }

  // conf.is3D() is assigned by the mol format parser based on the input
  // mol block designation, but also taking into account the presence of
  // non-null Z coordinates or stereobonds.
  //
  // the following test is in this sense probably redundant, but it's still
  // implemented in case molecules are built by other means.

  double max_absz {};
  for (const auto & p : conf.getPositions()) {
    max_absz = std::max(std::abs(p.z), max_absz);
  }

  if (max_absz > threshold) {
    errors.emplace_back("ERROR: [Is2DValidation] Molecule has non-null Z coordinates");
    if (!reportAllFailures) {
      return errors;
    }
  }

  if (conf.getNumAtoms() < 2) {
    // there is nothing else to check here, if there is at most one atom.
    return errors;
  }

  // verify that the atoms are not all in the same position (this often happens
  // because no coordinates were assigned and all atoms appear to be placed in
  // the origin)
  
  double min_x {}, max_x {};
  double min_y {}, max_y {};
  for (const auto & p : conf.getPositions()) {
    min_x = std::min(p.x, min_x);
    max_x = std::max(p.x, max_x);
    min_y = std::min(p.y, min_y);
    max_y = std::max(p.y, max_y);
  }
  auto delta_x = max_x - min_x;
  auto delta_y = max_y - min_y;
  auto max_delta = std::max(delta_x, delta_y);
  
  if (max_delta < threshold) {
    errors.emplace_back("ERROR: [Is2DValidation] All atoms have the same (x,y) coordinates");
    if (!reportAllFailures) {
      return errors;
    }
  }

  return errors;
}

namespace {
  double medianSquaredBondLength(const ROMol &mol, const Conformer &conf) {
    double median {};
    unsigned int numBonds = mol.getNumBonds();
    if (numBonds) {
      std::vector<double> values;
      values.reserve(numBonds);
      for (const auto & bond : mol.bonds()) {
        auto p1 = conf.getAtomPos(bond->getBeginAtomIdx());
        auto p2 = conf.getAtomPos(bond->getEndAtomIdx());
        auto dx = p2.x - p1.x;
        auto dy = p2.y - p1.y;
        auto value = dx*dx + dy*dy;
        values.push_back(value);
      }
      std::sort(values.begin(), values.end());
      if (numBonds % 2) {
        median = values[numBonds/2];
      }
      else {
        median = 0.5*(values[numBonds/2-1]+values[numBonds/2]);
      }
    }
    return median;
  }
}

std::vector<ValidationErrorInfo> Layout2DValidation::validate(
      const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  if (!mol.getNumConformers()) {
    errors.emplace_back("ERROR: [Layout2DValidation] Molecule has no conformer");
    return errors;
  }

  const auto & conf = mol.getConformer();
  unsigned int natoms = conf.getNumAtoms();

  if (natoms < 2) {
    // there is nothing to check here, if there is only one atom.
    return errors;
  }

  // compute threshold values for the squared atom-atom or atom-bond
  // distance and for the maximum bond length using the median squared
  // bond length as reference. 
  auto reference = medianSquaredBondLength(mol, conf);
  auto atomClashThreshold = clashLimit*clashLimit*reference;
  auto bondLengthThreshold = bondLengthLimit*bondLengthLimit*reference;

  for (unsigned int i = 0; i < natoms - 1; ++i) {
    const auto & pi = conf.getAtomPos(i);
    for (unsigned int j = i + 1; j < natoms; ++j) {
      const auto & pj = conf.getAtomPos(j);
      auto dx = pj.x - pi.x;
      auto dy = pj.y - pi.y;
      auto d2 = dx*dx + dy*dy;
      if (d2 < atomClashThreshold) {
        errors.push_back(
          "ERROR: [Layout2DValidation] atom " + std::to_string(i+1)
          + " too close to atom " + std::to_string(j+1));
        if (!reportAllFailures) {
          return errors;
        }
      }
    }
  }

  for (auto bond : mol.bonds()) {
    unsigned int i = bond->getBeginAtomIdx();
    const auto & pi = conf.getAtomPos(i);
    unsigned int j = bond->getEndAtomIdx();
    const auto & pj = conf.getAtomPos(j);

    auto ll = (pi.x - pj.x)*(pi.x - pj.x) + (pi.y - pj.y)*(pi.y - pj.y);
    if (ll > bondLengthThreshold) {
      errors.push_back(
        "ERROR: [Layout2DValidation] length of bond " + std::to_string(bond->getIdx()+1)
        + " between atoms " + std::to_string(i+1) + " and " + std::to_string(j+1)
        + " exceeds a configured limit");
      if (!reportAllFailures) {
        return errors;
      }
    }

    for (unsigned int k = 0; k < natoms; ++k) {
      if (k == i || k ==j) {
        continue;
      }
      const auto & pk = conf.getAtomPos(k);
      /*
               k
              /
            r/
            /
           /
          i---------------j
                  b
       */
      auto rr = (pk.x - pi.x)*(pk.x - pi.x) + (pk.y - pi.y)*(pk.y - pi.y);
      auto bb = (pj.x - pi.x)*(pj.x - pi.x) + (pj.y - pi.y)*(pj.y - pi.y);
      auto rb = (pk.x - pi.x)*(pj.x - pi.x) + (pk.y - pi.y)*(pj.y - pi.y);
      static constexpr double EPS {1.e-7}; // prevent dividing by zero in extreme cases
      auto kb = (rr*bb - rb*rb)/(bb+EPS);
      if (
        rb >= 0. &&   /* cos alpha > 0 */
        rb <= bb  &&  /* projection of r onto b does not exceed b */
        kb < atomClashThreshold /* distance from bond < limit */
        ) {
        errors.push_back(
          "ERROR: [Layout2DValidation] atom " + std::to_string(k+1)
          + " too close to bond " + std::to_string(bond->getIdx()+1));
        if (!reportAllFailures) {
          return errors;
        }
      }
    }
  }

  return errors;
}

namespace {
  bool isPotentialStereoCenter(const Atom * atom) {
    auto degree = atom->getDegree();
    auto atomicNum = atom->getAtomicNum();
    return (degree > 2 && degree <= 4 && (
        // the original set of atomic species recognized by STRUCHK
        // as stereogenic is C, N, O, Si, P, S
        atomicNum == 6 ||
        atomicNum == 7 ||
        // STRUCHK is documented to be based on the MDL guidelines, and
        // there is some evidence that oxygen was listed as a stereogenic species
        // in some old release of the MDL User Guide.
        // As of today, stereogenic oxygen atoms do apparently exist, but I'm not
        // really sure it is useful to include this possibility here.
        atomicNum == 8 ||
        atomicNum == 14 ||
        atomicNum == 15 ||
        atomicNum == 16
      ));
  }

  bool hasStereoBond(const ROMol &mol, const Atom * atom) {
    for (auto bond: mol.atomBonds(atom)) {
      if (atom != bond->getBeginAtom()) {
        continue;
      }
      auto bondDir = bond->getBondDir();
      if (
        bondDir == Bond::BondDir::BEGINDASH ||
        bondDir == Bond::BondDir::BEGINWEDGE ||
        bondDir == Bond::BondDir::UNKNOWN
        ) {
        return true;
      }
    }
    return false;
  }

  struct BondDirCount {
    unsigned int wedge = 0;
    unsigned int dash = 0;
    unsigned int unknown = 0;
    unsigned int other = 0;
  };

  struct NeighborsInfo {
    NeighborsInfo(const ROMol &mol, const Atom * atom);
    std::vector<const Bond *> bonds;
    BondDirCount dirCount;
  };

  NeighborsInfo::NeighborsInfo(const ROMol &mol, const Atom * atom)
  {
    for (auto bond: mol.atomBonds(atom)) {
      bonds.push_back(bond);
      Bond::BondDir dir = bond->getBondDir();
      switch (dir) {
      case Bond::BondDir::BEGINDASH:
        if (bond->getBeginAtom() == atom) {
          ++dirCount.dash;
        }
        break;
      case Bond::BondDir::BEGINWEDGE:
        if (bond->getBeginAtom() == atom) {
          ++dirCount.wedge;
        }
        break;
      case Bond::BondDir::UNKNOWN:
        if (bond->getBeginAtom() == atom) {
          ++dirCount.unknown;
        }
        break;
      case Bond::BondDir::NONE:
        // ok, ignore
        break;
      default:
        ++dirCount.other;
      }
    }

    const auto & conf = mol.getConformer();
    auto p = conf.getAtomPos(atom->getIdx());

    auto degree = bonds.size();
    std::vector<double> angles(degree);

    auto bond0 = bonds[0];
    auto atom0 = bond0->getOtherAtom(atom);
    auto v0 = conf.getAtomPos(atom0->getIdx()) - p;
    for (unsigned int n=1; n < degree; ++n) {
      auto bondn = bonds[n];
      auto atomn = bondn->getOtherAtom(atom);
      auto vn = conf.getAtomPos(atomn->getIdx()) - p;
      angles[n] = v0.signedAngleTo(vn);
    }

    // sort neighbors
    for (unsigned int i=2; i < degree; ++i) {
      for (unsigned int j=i; j>1; --j) {
        if (angles[j] < angles[j-1]) {
          auto bond = bonds[j];
          bonds[j] = bonds[j-1];
          bonds[j-1] = bond;
          auto angle = angles[j];
          angles[j] = angles[j-1];
          angles[j-1] = angle;
        }
        else {
          break;
        }
      }
    }
  }

  void check3CoordinatedStereo(
    const ROMol &mol, const Atom * atom,  const NeighborsInfo & neighborsInfo,
    bool /*reportAllFailures*/, std::vector<ValidationErrorInfo> & errors)
  {
    auto numStereoBonds = neighborsInfo.dirCount.dash + neighborsInfo.dirCount.wedge;

    if (numStereoBonds == 1) {
      // identify the stereo bond
      unsigned int i;
      for (i=0; i<3; ++i) {
        Bond::BondDir bondDir = neighborsInfo.bonds[i]->getBondDir();
        if (bondDir == Bond::BondDir::BEGINDASH || bondDir == Bond::BondDir::BEGINWEDGE) {
          break;
        }
      }
      // check for the colinearity of the stereocenter and the other two ligands.
      const auto & conf = mol.getConformer();
      auto p = conf.getAtomPos(atom->getIdx());
      auto atoma = neighborsInfo.bonds[(i+1)%3]->getOtherAtom(atom);
      auto va = conf.getAtomPos(atoma->getIdx()) - p;
      auto atomb = neighborsInfo.bonds[(i+2)%3]->getOtherAtom(atom);
      auto vb = conf.getAtomPos(atomb->getIdx()) - p;

      auto angle = va.angleTo(vb);

      static constexpr auto ANGLE_EPSILON = (M_PI*5./180.); // 5 degrees 
      if (angle < ANGLE_EPSILON || (M_PI - angle) < ANGLE_EPSILON) {
        errors.push_back(
          "ERROR: [StereoValidation] colinearity of non-stereo bonds at atom " + std::to_string(atom->getIdx()+1)
          );
      }
    }
    else {
      // configurations with multiple stereo bonds may be formally ambiguous or unambiguos
      // depending on their wedged/dashed direction and relative orientation on the plane.
      // those cases that are formally unambiguous are still most often discouraged or
      // also classified as not acceptable by IUPAC guidelines due to lack of clarity.

      // The AvalonTools' struchk implementation simply doesn't allow multiple stereo bonds
      // on stereo centers with 3 explicit ligands. The validations criteria for this sub-case
      // could be refined, but for now I'm keeping the same approach.
      errors.push_back(
        "ERROR: [StereoValidation] atom " + std::to_string(atom->getIdx()+1)
        + " has 3 explicit ligands and multiple stereo bonds"
        );
    }
  }

  void check4CoordinatedStereo(
    const ROMol &mol, const Atom * atom,  const NeighborsInfo & neighborsInfo,
    bool reportAllFailures, std::vector<ValidationErrorInfo> & errors)
  {
    if (neighborsInfo.dirCount.dash > 2 || neighborsInfo.dirCount.wedge > 2) {
      // this condition would anyway trigger an "adjacent bonds with like orientation"
      // alert, but this test could be clearer / more explicit.
      errors.push_back(
        "ERROR: [StereoValidation] atom " + std::to_string(atom->getIdx()+1)
        + " has too many stereo bonds with like orientation"
        );
      if (!reportAllFailures) {
        return;
      }
    }

    for (unsigned int i=0; i<2; ++i) {
      if ((neighborsInfo.bonds[i]->getBondDir() == Bond::BondDir::BEGINDASH &&
           neighborsInfo.bonds[i+2]->getBondDir() == Bond::BondDir::BEGINWEDGE)
          ||
          (neighborsInfo.bonds[i]->getBondDir() == Bond::BondDir::BEGINWEDGE &&
           neighborsInfo.bonds[i+2]->getBondDir() == Bond::BondDir::BEGINDASH)) {
        errors.push_back(
          "ERROR: [StereoValidation] atom " + std::to_string(atom->getIdx()+1)
          + " has opposing stereo bonds with different up/down orientation"
          );
        if (!reportAllFailures) {
          return;
        }
      }    
    }

    for (unsigned int i=0; i<4; ++i) {
      if ((neighborsInfo.bonds[i]->getBondDir() == Bond::BondDir::BEGINDASH &&
           neighborsInfo.bonds[(i+1)%4]->getBondDir() == Bond::BondDir::BEGINDASH)
          ||
          (neighborsInfo.bonds[i]->getBondDir() == Bond::BondDir::BEGINWEDGE &&
           neighborsInfo.bonds[(i+1)%4]->getBondDir() == Bond::BondDir::BEGINWEDGE)) {
        errors.push_back(
          "ERROR: [StereoValidation] atom " + std::to_string(atom->getIdx()+1)
          + " has adjacent stereo bonds with like orientation"
          );
        if (!reportAllFailures) {
          return;
        }
        // it doesn't make sense to output this alert multiple times for the same atom
        // we therefore exit the loop also when reportAllFailures is not set.
        break; 
      }    
    }

    if (neighborsInfo.dirCount.dash + neighborsInfo.dirCount.wedge == 1) {

      // there is only one wedged/dashed bond. check for 'umbrellas' and
      // other geometric violations. we need the conformation here.
      const auto & conf = mol.getConformer();

      // identify the bond index for the stereo bond with specified direction.
      for (unsigned int i=0; i<4; ++i) {
        Bond::BondDir bondDir = neighborsInfo.bonds[i]->getBondDir();
        if (bondDir == Bond::BondDir::BEGINDASH || bondDir == Bond::BondDir::BEGINWEDGE) {
          // count how many of the other bonds lie on the opposite half-plane, i.e.
          // form an angle > pi/4 with the stereo bond (compute the dot-product of the
          // corresponding vectors, and check if it's negative).
          unsigned int opposed = 0;
          auto p = conf.getAtomPos(atom->getIdx());
          auto bondi = neighborsInfo.bonds[i];
          auto atomi = bondi->getOtherAtom(atom);
          auto vi = conf.getAtomPos(atomi->getIdx()) - p;
          for (unsigned int j=0; j<4; ++j) {
            if (j==i) {
              continue;
            }
            auto bondj = neighborsInfo.bonds[j];
            auto atomj = bondj->getOtherAtom(atom);
            auto vj = conf.getAtomPos(atomj->getIdx()) - p;
            if (vi.x*vj.x + vi.y*vj.y < 0.) {
              ++opposed;
            }
          }
          if (opposed == 3) {
            errors.push_back(
              "ERROR: [StereoValidation] atom " + std::to_string(atom->getIdx()+1)
              + " has a potentially ambiguous representation: all non-stereo bonds"
              + " opposite to the only stereo bond"
              );
          }
          if (!reportAllFailures) {
            return;
          }
          // there is only one stereo bond, which means we can exit the
          // outer loop on the first execution of this block.
          break;
        }
      }

      // check for collinearity violations and/or cases where the
      // the middle non-stereo bond is badly positioned (i.e., too short
      // compared to the other two on its sides).
      for (unsigned int i=0; i<4; i++) {
        Bond::BondDir bondDir = neighborsInfo.bonds[i]->getBondDir();
        if (bondDir == Bond::BondDir::BEGINDASH || bondDir == Bond::BondDir::BEGINWEDGE) {
          auto j = (i+1) % 4;
          auto k = (i+2) % 4;
          auto l = (i+3) % 4;
          auto atomj = neighborsInfo.bonds[j]->getOtherAtom(atom);
          auto atomk = neighborsInfo.bonds[k]->getOtherAtom(atom);
          auto atoml = neighborsInfo.bonds[l]->getOtherAtom(atom);
          auto pj = conf.getAtomPos(atomj->getIdx());
          auto pk = conf.getAtomPos(atomk->getIdx());
          auto pl = conf.getAtomPos(atoml->getIdx());
          auto v1 = pj - pk;
          auto v2 = pl - pk;
          auto angle = v1.signedAngleTo(v2);
          if (angle < 185.*M_PI/180.) {
            errors.push_back(
              "ERROR: [StereoValidation] colinearity or triangle rule violation of "
              "non-stereo bonds at atom " + std::to_string(atom->getIdx()+1) /* +
              " due to angle formed by ("  +
              std::to_string(atomj->getIdx()+1) + "," +
              std::to_string(atomk->getIdx()+1) + "," +
              std::to_string(atoml->getIdx()+1) + ")" */
            );
            if (!reportAllFailures) {
              return;
            }
          }
          // there is only one stereo bond, which means we can exit the
          // outer loop on the first execution of this block.
          break;
        }
      }
    }

  }

  void checkStereo(
    const ROMol &mol, const Atom * atom,  bool reportAllFailures,
    std::vector<ValidationErrorInfo> & errors)
  {
    bool multipleBondFound {}, possibleAllene {};
    for (auto bond: mol.atomBonds(atom)) {
      auto bondType = bond->getBondType();
      if (bondType != Bond::BondType::SINGLE) {
        multipleBondFound = true;
        const Atom * otherAtom = bond->getOtherAtom(atom);
        if (otherAtom->getDegree() == 2) {
          int doubleBondCount {};
          for (auto otherBond: mol.atomBonds(otherAtom)) {
            if (otherBond->getBondType() == Bond::BondType::DOUBLE) {
              ++doubleBondCount;
            }
          }
          if (doubleBondCount == 2) {
            possibleAllene = true;
          }
        }
      }
    }
    auto atomicNum = atom->getAtomicNum();
    if (
      multipleBondFound && !possibleAllene && atomicNum != 15 && atomicNum != 16
      ) {
        errors.push_back(
          "ERROR: [StereoValidation] unexpected stereo bond found at unsaturated atom "
          + std::to_string(atom->getIdx()+1));
        return;
    }
    else if (multipleBondFound && atomicNum != 16) {
      // Other cases of unsaturated atoms (allenes, P compounds)
      // are not further validated at this time.
      return;
    }

    NeighborsInfo neighborsInfo(mol, atom);

    if (neighborsInfo.dirCount.other) {
      errors.push_back(
        "ERROR: [StereoValidation] one or more bonds incident to atom "
        + std::to_string(atom->getIdx()+1)
        + "have invalid direction settings");
      // this is an unlikely condition and it would make little sense to
      // continue the analysis also when reportAllFailures were set.
      return;
    }

    if (neighborsInfo.dirCount.unknown) {
      if (neighborsInfo.dirCount.dash || neighborsInfo.dirCount.wedge) {
        errors.push_back(
          "ERROR: [StereoValidation] atom " + std::to_string(atom->getIdx()+1)
          + "has both unknown and wedged/dashed stereo bonds.");
      }
      // else: if the only stereo bonds have either/unknown direction,
      // we can anyway return here.
      return;
    }

    if (atom->getDegree() == 3) {
      check3CoordinatedStereo(mol, atom, neighborsInfo, reportAllFailures, errors);
    }
    else {
      check4CoordinatedStereo(mol, atom, neighborsInfo, reportAllFailures, errors);
    }
  }
}

std::vector<ValidationErrorInfo> StereoValidation::validate(
      const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;

  for (auto atom: mol.atoms()) {
    bool stereoBondFound = hasStereoBond(mol, atom);
    if (stereoBondFound && isPotentialStereoCenter(atom)) {
      checkStereo(mol, atom, reportAllFailures, errors);
    }
    /* do not disallow stereo bonds on non-stereogenic centers because
     * they may have been used to describe some kind of atropisomerism.
     * the else branch below is for this reason commented-out.
     * note: no attempt is made to verify that this usage of stereo bonds
     * is well-defined/unambiguous.
     */
    /* else if (stereoBondFound) {
      errors.push_back(
        "ERROR: [StereoValidation] atom " + std::to_string(atom->getIdx()+1)
        + " has stereo bonds, but it doesn't seem to be a stereogenic center");
    }*/
    if (!errors.empty() && !reportAllFailures) {
      break;
    }
  }

  return errors;
}

std::vector<ValidationErrorInfo> validateSmiles(const std::string &smiles) {
  RWMOL_SPTR mol(SmilesToMol(smiles));
  if (!mol) {
    std::string message =
        "SMILES Parse Error: syntax error for input: " + smiles;
    throw ValueErrorException(message);
  }

  MolVSValidation vm;
  std::vector<ValidationErrorInfo> errors = vm.validate(*mol, true);

  return errors;
}

}  // namespace MolStandardize
}  // namespace RDKit
