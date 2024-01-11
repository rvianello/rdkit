//
//  Copyright (C) 2023 Novartis Biomedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Pipeline.h"
#include "Validate.h"
#include "Metal.h"
#include "Normalize.h"
#include "Charge.h"
#include "Fragment.h"
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Chirality.h>

namespace RDKit {
namespace MolStandardize {

void PipelineResult::append(PipelineStatus newStatus, const std::string & info)
{
  status = static_cast<PipelineStatus>(status | newStatus);
  log.push_back({newStatus, info});
}

PipelineResult Pipeline::run(const std::string & molblock)
{
  PipelineResult result;
  result.status = NO_ERROR;
  result.stage = STARTED;
  result.inputMolBlock = molblock;
  result.outputMolBlock = molblock;

  // parse the molblock into an RWMol instance
  std::unique_ptr<RWMol> mol = parse(molblock, result);
  if (!mol || (result.status != NO_ERROR && !options.reportAllFailures)) {
    return result;
  }

  // validate the structure
  validate(*mol, result);
  if (result.status != NO_ERROR && !options.reportAllFailures) {
    return result;
  }

  // standardize/normalize
  mol = standardize(*mol, result);
  if (!mol || (result.status != NO_ERROR && !options.reportAllFailures)) {
    return result;
  }

  // serialize as MolBlock
  serialize(*mol, result);
  if (result.status != NO_ERROR && !options.reportAllFailures) {
    return result;
  }

  result.stage = COMPLETED;

  return result;
}

std::unique_ptr<RWMol> Pipeline::parse(const std::string & molblock, PipelineResult & result)
{
  result.stage = PARSING_INPUT;

  // we don't want to sanitize the molecule at this stage
  static constexpr bool sanitize {false};
  // Hs are anyways not removed if the mol is not sanitized
  static constexpr bool removeHs {false};

  // strict parsing is configurable via the pipeline options
  const bool strictParsing {options.strictParsing};

  std::unique_ptr<RWMol> mol {};

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

void Pipeline::validate(const ROMol & mol, PipelineResult & result)
{
  result.stage = VALIDATION;

  auto applyValidation = [&mol, &result, this](const ValidationMethod & v, PipelineStatus status) -> bool {
    auto errors = v.validate(mol, options.reportAllFailures);
    for (const auto & error : errors) {
      result.append(status, error.what());
    }
    return errors.empty();
  };

  // check the number of atoms and valence status
  RDKitValidation rdkitValidation;
  if (!applyValidation(rdkitValidation, BASIC_VALIDATION_ERROR) && !options.reportAllFailures) {
    return;
  }

  // check for undesired features in the input molecule
  // TBD

  // verify that the input is a 2D structure
  Is2DValidation is2DValidation(options.is2DZeroThreshold);
  if (!applyValidation(is2DValidation, IS2D_VALIDATION_ERROR) && !options.reportAllFailures) {
    return;
  }

  // validate the 2D layout (check for clashing atoms and abnormally long bonds)
  Layout2DValidation layout2DValidation(options.atomClashLimit, options.bondLengthLimit);
  if (!applyValidation(layout2DValidation, LAYOUT2D_VALIDATION_ERROR) && !options.reportAllFailures) {
    return;
  }

  // verify that the specified stereochemistry is formally correct
  StereoValidation stereoValidation;
  if (!applyValidation(stereoValidation, STEREO_VALIDATION_ERROR) && !options.reportAllFailures) {
    return;
  }
}

std::unique_ptr<RWMol> Pipeline::standardize(RWMol & original, PipelineResult & result)
{
  result.stage = STANDARDIZATION;

  std::unique_ptr<RWMol> mol {new RWMol(original)};

  // Prepare the mol for standardization
  try {
    constexpr unsigned int sanitizeOps = MolOps::SANITIZE_ALL
      ^ MolOps::SANITIZE_CLEANUP
      ^ MolOps::SANITIZE_CLEANUP_ORGANOMETALLICS
      ^ MolOps::SANITIZE_PROPERTIES
      ;
    unsigned int failedOp = 0;
    MolOps::sanitizeMol(*mol, failedOp, sanitizeOps);
  }
  catch (MolSanitizeException &) {
    result.append(
      PREPARE_STANDARDIZATION_ERROR,
      "ERROR: [Standardization] An unexpected error occurred while preparing the molecule for standardization.");
    return std::unique_ptr<RWMol> {};
  }

  // bonding to metals
  try {
    MetalDisconnector metalDisconnector;
    metalDisconnector.disconnectInPlace(*mol);
  }
  catch (...) {
    result.append(
      METAL_STANDARDIZATION_ERROR,
      "ERROR: [Standardization] An unexpected error occurred while processing the bonding of metal species.");
    return std::unique_ptr<RWMol> {};
  }

  // functional groups
  try {
    Normalizer normalizer;
    normalizer.normalizeInPlace(*mol);
  }
  catch (...) {
    result.append(
      NORMALIZER_STANDARDIZATION_ERROR, 
      "ERROR: [Standardization] An unexpected error occurred while normalizing the representation of some functional groups");
    return std::unique_ptr<RWMol> {};
  }

  // rearrange the protonation status
  // TODO: is this useful/necessary, considering the Uncharger applied at the end?
  try {
    Reionizer reionizer;
    reionizer.reionizeInPlace(*mol);
  }
  catch (...) {
    result.append(
      PROTONATION_STANDARDIZATION_ERROR, 
      "ERROR: [Standardization] An unexpected error occurred while reassigning the formal charges");
    return std::unique_ptr<RWMol> {};
  }

  // keep the largest fragment
  try {
    LargestFragmentChooser fragmentChooser;
    std::unique_ptr<ROMol> largestFragment {fragmentChooser.choose(*mol)};
    mol.reset(static_cast<RWMol *>(largestFragment.release()));
  }
  catch (...) {
    result.append(
      FRAGMENT_STANDARDIZATION_ERROR, 
      "ERROR: [Standardization] An unexpected error occurred while removing the disconnected fragments");
  }

  // overall charge status
  try {
    Uncharger uncharger;
    uncharger.unchargeInPlace(*mol);
  }
  catch (...) {
    result.append(
      CHARGE_STANDARDIZATION_ERROR, 
      "ERROR: [Standardization] An unexpected error occurred while normalizing the total charge status");
    return std::unique_ptr<RWMol> {};
  }

  // updating the property cache was observed to be required, in order to clear the explicit valence
  // property (CTab's VAL) from deprotonated quaternary nitrogens, that would otherwise persist in the
  // output MolBlock and would result in an invalid molecule.
  //
  // needsUpdatePropertyCache() returns false 
  //if (mol->needsUpdatePropertyCache()) {
  if (true) {
    mol->updatePropertyCache(false);
  }

  return mol;
}

void Pipeline::serialize(const ROMol & mol, PipelineResult & result)
{
  result.stage = SERIALIZING_OUTPUT;

  try {
    if (!options.ouputV2000) {
      result.outputMolBlock = MolToV3KMolBlock(mol);
    }
    else if (mol.getNumAtoms() > 999 || mol.getNumBonds() > 999) {
      result.append(OUTPUT_ERROR, "ERROR: [Output] Molecule is too large for V2000 format.");
    }
    else {
      result.outputMolBlock = MolToMolBlock(mol);
    }
  }
  catch (const std::exception & e) {
    result.append(
      OUTPUT_ERROR,
      "ERROR: [Output] Can't write molecule to output format: " + std::string(e.what()));
  }
  catch (...) {
    result.append(
      OUTPUT_ERROR,
      "ERROR: [Output] An unexpected error occurred.");
  }
}

}
}
