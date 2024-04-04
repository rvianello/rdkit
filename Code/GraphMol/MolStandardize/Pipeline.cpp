//
//  Copyright (C) 2023 Novartis Biomedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <cmath>
#include <regex>
#include <sstream>
#include "Pipeline.h"
#include "Validate.h"
#include "Metal.h"
#include "Normalize.h"
#include "Charge.h"
#include "Fragment.h"
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Chirality.h>

namespace RDKit {
namespace MolStandardize {

void PipelineResult::append(PipelineStatus newStatus, const std::string & info)
{
  status = static_cast<PipelineStatus>(status | newStatus);
  log.push_back({newStatus, info});
}

PipelineResult Pipeline::run(const std::string & molblock) const
{
  PipelineResult result;
  result.status = NO_EVENT;
  result.inputMolBlock = molblock;

  // parse the molblock into an RWMol instance
  result.stage = PARSING_INPUT;
  RWMOL_SPTR mol = parse(molblock, result);
  if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures)) {
    return result;
  }

  RWMOL_SPTR_PAIR output;

  if (mol->getNumAtoms() == 0 && options.allowEmptyMolecules) {
    output = {mol, mol};
  }
  else {
    // input sanitization + cleanup
    result.stage = PREPARE_FOR_VALIDATION;
    mol = prepareForValidation(mol, result);
    if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures)) {
      return result;
    }

    // validate the structure
    result.stage = VALIDATION;
    mol = validate(mol, result);
    if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures)) {
      return result;
    }

    // re-read and sanitize the validated structure
    result.stage = PREPARE_FOR_STANDARDIZATION;
    mol = parse(molblock, result);
    if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures)) {
      return result;
    }
    mol = prepareForStandardization(mol, result);
    if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures)) {
      return result;
    }

    // standardize/normalize
    result.stage = STANDARDIZATION;
    mol = standardize(mol, result);
    if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures)) {
      return result;
    }
    output = makeParent(mol, result);
    if (!output.first || !output.second
        || ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures)) {
      return result;
    }
  }

  // serialize as MolBlocks
  result.stage = SERIALIZING_OUTPUT;
  serialize(output, result);
  if ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures) {
    return result;
  }

  result.stage = COMPLETED;

  return result;
}

RWMOL_SPTR Pipeline::parse(const std::string & molblock, PipelineResult & result) const
{
  // we don't want to sanitize the molecule at this stage
  static constexpr bool sanitize {false};
  // Hs wouldn't be anyway removed if the mol is not sanitized
  static constexpr bool removeHs {false};

  // strict parsing is configurable via the pipeline options
  const bool strictParsing {options.strictParsing};

  RWMOL_SPTR mol {};

  try {
    mol.reset(MolBlockToMol(molblock, sanitize, removeHs, strictParsing));
  }
  catch (FileParseException & e) {
    result.append(INPUT_ERROR, e.what());
  }

  if (!mol) {
    result.append(
      INPUT_ERROR,
      "Could not instantiate a valid molecule from input"
    );
  }

  return mol;
}

RWMOL_SPTR Pipeline::prepareForValidation(RWMOL_SPTR mol, PipelineResult & result) const
{
  // Prepare the mol for validation.

  try {
    // convert to smiles and later check if the structure was modified
    auto reference = MolToSmiles(*mol);

    // The general intention is about validating the original input, and therefore
    // limit the sanitization to the minimum, but it's not very useful to record a
    // valence validation error for issues like a badly drawn nitro group that would
    // be later fixed during by the normalization step.
    //
    // Some sanitization also needs to be performed in order to assign the stereochemistry
    // (which needs to happen prior to reapplying the wedging, see below), and we need
    // to find radicals, in order to support the corresponding validation criterion.
    constexpr unsigned int sanitizeOps = (
      MolOps::SANITIZE_CLEANUP
      | MolOps::SANITIZE_SYMMRINGS
      | MolOps::SANITIZE_CLEANUP_ORGANOMETALLICS
      | MolOps::SANITIZE_FINDRADICALS
    );
    unsigned int failedOp = 0;
    MolOps::sanitizeMol(*mol, failedOp, sanitizeOps);

    auto smiles = MolToSmiles(*mol);
    if (reference != smiles) {
      result.append(
        SANITIZATION_APPLIED,
        "Some traits in the representation of the chemical structure were updated in a pre-validation cleanup step.");
    }

    // We want to restore the original MolBlock wedging, but this step may in some cases overwrite the
    // ENDDOWNRIGHT/ENDUPRIGHT info that describes the configuration of double bonds adjacent to stereocenters.
    // We therefore need to assign the double bond stereochemistry now, so that they are not marked as
    // EITHERDOUBLE in the output MolBlocks.
    MolOps::assignStereochemistry(*mol, true, true, true);
    Chirality::reapplyMolBlockWedging(*mol);
  }
  catch (MolSanitizeException &) {
    result.append(
      PREPARE_FOR_VALIDATION_ERROR,
      "An unexpected error occurred while preparing the molecule for validation.");
  }

  return mol;
}

namespace {
  // The error messages from the ValidationMethod classes include some metadata
  // in a string prefix that are not particularly useful within the context of this
  // Pipeline. The function below removes that prefix.
  static const std::regex prefix("^(ERROR|INFO): \\[.+\\] ");
  std::string removeErrorPrefix(const std::string & message) {
    return std::regex_replace(message, prefix, "");
  }
}

RWMOL_SPTR Pipeline::validate(RWMOL_SPTR mol, PipelineResult & result) const
{
  auto applyValidation = [&mol, &result, this](const ValidationMethod & v, PipelineStatus status) -> bool {
    auto errors = v.validate(*mol, options.reportAllFailures);
    for (const auto & error : errors) {
      result.append(status, removeErrorPrefix(error));
    }
    return errors.empty();
  };

  // check for undesired features in the input molecule (e.g., query atoms/bonds)
  FeaturesValidation featuresValidation(options.allowEnhancedStereo, options.allowAromaticBondType);
  if (!applyValidation(featuresValidation, FEATURES_VALIDATION_ERROR) && !options.reportAllFailures) {
    return mol;
  }

  // check the number of atoms and valence status
  RDKitValidation rdkitValidation;
  if (!applyValidation(rdkitValidation, BASIC_VALIDATION_ERROR) && !options.reportAllFailures) {
    return mol;
  }

  // disallow radicals
  DisallowedRadicalValidation radicalValidation;
  if (!applyValidation(radicalValidation, BASIC_VALIDATION_ERROR) && !options.reportAllFailures) {
    return mol;
  }

  // validate the isotopic numbers (if any are specified)
  IsotopeValidation isotopeValidation(true);
  if (!applyValidation(isotopeValidation, BASIC_VALIDATION_ERROR) && !options.reportAllFailures) {
    return mol;
  }

  // verify that the input is a 2D structure
  Is2DValidation is2DValidation(options.is2DZeroThreshold);
  if (!applyValidation(is2DValidation, IS2D_VALIDATION_ERROR) && !options.reportAllFailures) {
    return mol;
  }

  // validate the 2D layout (check for clashing atoms and abnormally long bonds)
  Layout2DValidation layout2DValidation(
    options.atomClashLimit, options.bondLengthLimit,
    options.allowLongBondsInRings, options.allowAtomBondClashExemption,
    options.minMedianBondLength);
  if (!applyValidation(layout2DValidation, LAYOUT2D_VALIDATION_ERROR) && !options.reportAllFailures) {
    return mol;
  }

  // verify that the specified stereochemistry is formally correct
  StereoValidation stereoValidation;
  if (!applyValidation(stereoValidation, STEREO_VALIDATION_ERROR) && !options.reportAllFailures) {
    return mol;
  }

  return mol;
}

RWMOL_SPTR Pipeline::prepareForStandardization(RWMOL_SPTR mol, PipelineResult & result) const
{
  // Prepare the mol for standardization.

  try {
    MolOps::sanitizeMol(*mol);
  }
  catch (MolSanitizeException &) {
    result.append(
      PREPARE_FOR_STANDARDIZATION_ERROR,
      "An unexpected error occurred while preparing the molecule for standardization.");
  }

  return mol;
}

RWMOL_SPTR Pipeline::standardize(RWMOL_SPTR mol, PipelineResult & result) const
{
  auto smiles = MolToSmiles(*mol);
  auto reference = smiles;

  // bonding to metals
  try {
    MetalDisconnectorOptions mdOpts;
    mdOpts.allowPartialDisconnections = options.allowPartialDisconnections;
    MetalDisconnector metalDisconnector(mdOpts);
    std::unique_ptr<ROMol> metalNof {SmartsToMol(options.metalNof)};
    metalDisconnector.setMetalNof(*metalNof);
    std::unique_ptr<ROMol> metalNon {SmartsToMol(options.metalNon)};
    metalDisconnector.setMetalNon(*metalNon);
    metalDisconnector.disconnectInPlace(*mol);
  }
  catch (...) {
    result.append(
      METAL_STANDARDIZATION_ERROR,
      "An unexpected error occurred while processing the bonding of metal species.");
    return mol;
  }

  smiles = MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(METALS_DISCONNECTED, "One or more metal atoms were disconnected.");
  }
  reference = smiles;

  // functional groups
  try {
    std::unique_ptr<Normalizer> normalizer {};
    if (options.normalizerData.empty()) {
      normalizer.reset(new Normalizer);
    }
    else {
      std::istringstream sstr(options.normalizerData);
      normalizer.reset(new Normalizer(sstr, options.normalizerMaxRestarts));
    }
    normalizer->normalizeInPlace(*mol);
  }
  catch (...) {
    result.append(
      NORMALIZER_STANDARDIZATION_ERROR, 
      "An unexpected error occurred while normalizing the representation of some functional groups");
    return mol;
  }

  smiles = MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(
      NORMALIZATION_APPLIED,
      "The representation of some functional groups was adjusted.");
  }
  reference = smiles;

  // keep the largest fragment
  try {
    LargestFragmentChooser fragmentChooser;
    fragmentChooser.chooseInPlace(*mol);
  }
  catch (...) {
    result.append(
      FRAGMENT_STANDARDIZATION_ERROR, 
      "An unexpected error occurred while removing the disconnected fragments");
    return mol;
  }

  smiles = MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(
      FRAGMENTS_REMOVED,
      "One or more disconnected fragments (e.g., counterions) were removed.");
  }

  // The stereochemistry is not assigned until after we are done modifying the
  // molecular graph:
  MolOps::assignStereochemistry(*mol, true, true, true);
  // restore the original wedging from the MolBlock
  Chirality::reapplyMolBlockWedging(*mol);

  // scale the atoms coordinates
  // and make sure that z coords are set to 0 (some z coords may be non-null
  // albeit smaller than the validation threshold - these noisy coords may in some cases
  // also interfere with the perception of stereochemistry by some tools e.g., inchi)
  if (options.scaledMedianBondLength > 0. && mol->getNumConformers()) {
    auto & conf = mol->getConformer();
    double medianBondLength = sqrt(Layout2DValidation::squaredMedianBondLength(*mol, conf));
    if (medianBondLength > options.minMedianBondLength) {
      double scaleFactor = options.scaledMedianBondLength/medianBondLength;
      unsigned int natoms = conf.getNumAtoms();
      for (unsigned int i = 0; i < natoms; ++i) {
        auto pos = conf.getAtomPos(i)*scaleFactor;
        pos.z = 0.;
        conf.setAtomPos(i, pos);
      }
    }
  }


  return mol;
}

namespace {

  void removeHsAtProtonatedSites(RWMOL_SPTR mol) {
    boost::dynamic_bitset<> protons{mol->getNumAtoms(), 0};
    for (auto atom : mol->atoms()) {
      if (atom->getAtomicNum() != 1 || atom->getDegree() != 1) {
        continue;
      }
      for (auto neighbor : mol->atomNeighbors(atom)) {
        if (neighbor->getFormalCharge() > 0) {
          protons.set(atom->getIdx());
        }
      }
    }
    if (protons.any()) {
      for (int idx = mol->getNumAtoms() - 1; idx >= 0; --idx) {
        if (!protons[idx]) {
          continue;
        }
        auto atom = mol->getAtomWithIdx(idx);
        for (auto bond : mol->atomBonds(atom)) {
          auto neighbor = bond->getOtherAtom(atom);
          neighbor->setNumExplicitHs(neighbor->getNumExplicitHs() + 1);
          break; // there are no other bonds anyways
        }
        mol->removeAtom(atom);
      }
      mol->updatePropertyCache(false);
    }
  }
}

Pipeline::RWMOL_SPTR_PAIR Pipeline::makeParent(RWMOL_SPTR mol, PipelineResult & result) const
{
  auto reference = MolToSmiles(*mol);

  RWMOL_SPTR parent {new RWMol(*mol)};

  // overall charge status
  try {
    // The Uncharger implementation wouldn't identify the positively
    // charged sites with adjacent explicit Hs correctly (it's a quite
    // unlikely configuration, but potentially possible considering that
    // the pipeline operates on unsanitized input).
    //
    // If present, these Hs are therefore removed from the molecular graph
    // prior to neutralization.
    removeHsAtProtonatedSites(parent);

    static const bool canonicalOrdering = false;
    static const bool force = true;
    Uncharger uncharger(canonicalOrdering, force);
    uncharger.unchargeInPlace(*parent);
  }
  catch (...) {
    result.append(
      CHARGE_STANDARDIZATION_ERROR, 
      "An unexpected error occurred while normalizing the compound's charge status");
    return {{}, {}};
  }

  // Check if `mol` was submitted in a suitable ionization state
  int parentCharge {};
  for (auto atom: parent->atoms()) {
    parentCharge += atom->getFormalCharge();
  }

  int molCharge {};
  for (auto atom: mol->atoms()) {
    molCharge += atom->getFormalCharge();
  }

  // If mol is neutral or in a protonation state that partially or fully
  // balances the non-neutralizable charged sites in the parent structure,
  // then mol is accepted. Otherwise, it is replaced by its parent.
  if ((molCharge > 0 && molCharge > parentCharge) ||
      (molCharge < 0 && molCharge < parentCharge)) {
    mol = parent;
  }

  auto smiles = MolToSmiles(*mol);
  if (smiles != reference) {
    result.append(PROTONATION_CHANGED, "The protonation state was adjusted.");
  }
  reference = smiles;

  // updating the property cache was observed to be required, in order to clear the explicit valence
  // property (CTab's VAL) from deprotonated quaternary nitrogens, that could be specified in the
  // original input and otherwise persisted in the output MolBlock, resulting in an invalid molecule.
  //
  // needsUpdatePropertyCache() returns false 
  //if (mol->needsUpdatePropertyCache()) {
  mol->updatePropertyCache(false);
  parent->updatePropertyCache(false);

  return {mol, parent};
}

void Pipeline::serialize(RWMOL_SPTR_PAIR output, PipelineResult & result) const
{
  const ROMol & outputMol = *output.first;
  const ROMol & parentMol = *output.second;
  
  try {
    if (!options.ouputV2000) {
      result.outputMolBlock = MolToV3KMolBlock(outputMol);
      result.parentMolBlock = MolToV3KMolBlock(parentMol);
    }
    else if (
        outputMol.getNumAtoms() > 999 || outputMol.getNumBonds() > 999 ||
        parentMol.getNumAtoms() > 999 || parentMol.getNumBonds() > 999
        ) {
      result.append(OUTPUT_ERROR, "The molecule is too large for V2000 format.");
    }
    else {
      result.outputMolBlock = MolToMolBlock(outputMol);
      result.parentMolBlock = MolToMolBlock(parentMol);
    }
  }
  catch (const std::exception & e) {
    result.append(
      OUTPUT_ERROR,
      "Can't write molecule to output format: " + std::string(e.what()));
  }
  catch (...) {
    result.append(
      OUTPUT_ERROR,
      "An unexpected error occurred while serializing the output structures.");
  }
}

}
}
