//
//  Copyright (C) 2023 Novartis Biomedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

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
    result.stage = SANITIZATION;
    mol = sanitize(mol, result);
    if (!mol || ((result.status & PIPELINE_ERROR) != NO_EVENT && !options.reportAllFailures)) {
      return result;
    }

    // validate the structure
    result.stage = VALIDATION;
    mol = validate(mol, result);
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
  else {
    Chirality::reapplyMolBlockWedging(*mol);
  }

  return mol;
}

RWMOL_SPTR Pipeline::sanitize(RWMOL_SPTR mol, PipelineResult & result) const
{
  // Prepare the mol for validation and standardization.
  //
  // The general intention is about validating the original input, and therefore
  // limit the sanitization to the minimum, but it's not very useful to record a
  // valence validation error for issues like a badly drawn nitro group that will
  // be later fixed during by the normalization step. Moreover, some standardization
  // components do not operate correctly if sanitization is not performed at all.

  try {
    // convert to smiles and later check if the structure was modified
    auto reference = MolToSmiles(*mol);
    // TODO: consider extending this pre-validation step to other
    // sanitization flags that may be needed/useful, but do not introduce
    // significant changes in the input.
    constexpr unsigned int sanitizeOps = MolOps::SANITIZE_CLEANUP;
    unsigned int failedOp = 0;
    MolOps::sanitizeMol(*mol, failedOp, sanitizeOps);
    auto smiles = MolToSmiles(*mol);
    if (reference != smiles) {
      result.append(
        SANITIZATION_APPLIED,
        "The representation of some functional groups was modified in a pre-validation cleanup step.");
    }
  }
  catch (MolSanitizeException &) {
    result.append(
      SANITIZATION_ERROR,
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
  FeaturesValidation featuresValidation(options.allowEnhancedStereo);
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
  Layout2DValidation layout2DValidation(options.atomClashLimit, options.bondLengthLimit);
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

RWMOL_SPTR Pipeline::standardize(RWMOL_SPTR mol, PipelineResult & result) const
{
  auto smiles = MolToSmiles(*mol);
  auto reference = smiles;

  // bonding to metals
  try {
    MetalDisconnector metalDisconnector;
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

  return mol;
}

Pipeline::RWMOL_SPTR_PAIR Pipeline::makeParent(RWMOL_SPTR mol, PipelineResult & result) const
{
  auto reference = MolToSmiles(*mol);

  RWMOL_SPTR parent {new RWMol(*mol)};

  // overall charge status
  try {
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

  int parentCharge {};
  for (auto atom: parent->atoms()) {
    parentCharge += atom->getFormalCharge();
  }

  if (parentCharge <0) {
    // this is actually unexpected
    result.append(
      CHARGE_STANDARDIZATION_ERROR,
      "Could not produce a valid uncharged structure");
    return {{}, {}};
  }

  // Check if `mol` was submitted in a suitable ionization state
  int molCharge {};
  for (auto atom: mol->atoms()) {
    molCharge += atom->getFormalCharge();
  }

  // If mol is neutral or in a protonation state that partially or fully
  // balances the non-neutralizable positively charged sites in the parent
  // structure, then mol is accepted. Otherwise, it is replaced by its parent.
  if (molCharge > parentCharge || molCharge < 0) {
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
