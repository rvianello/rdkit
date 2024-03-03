//
//  Copyright (C) 2023 Novartis Biomedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_MOLSTANDARDIZE_PIPELINE_H
#define RD_MOLSTANDARDIZE_PIPELINE_H
#include <RDGeneral/export.h>
#include <GraphMol/RWMol.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace RDKit {

namespace MolStandardize {

struct RDKIT_MOLSTANDARDIZE_EXPORT PipelineOptions {
  // parsing
  bool strictParsing {false};

  // validation
  bool reportAllFailures {true};
  bool allowEmptyMolecules {false};
  bool allowEnhancedStereo {false};
  bool allowAromaticBondType {false};
  double is2DZeroThreshold {1e-3};
  double atomClashLimit {0.03};
  double bondLengthLimit {50.};
  bool allowLongBondsInRings {true};

  // cleanup/standardization
  // metal disconnector options
  std::string metalNof {"[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra]~[#7,#8,F]"};
  std::string metalNon {};
  // normalizer options
  std::string normalizerData {};
  unsigned int normalizerMaxRestarts {200};

  // serialization
  bool ouputV2000 {false};
};

enum RDKIT_MOLSTANDARDIZE_EXPORT PipelineStatus {
  NO_EVENT=0,
  INPUT_ERROR=(1<<0),
  PREPARE_FOR_VALIDATION_ERROR=(1<<1),
  FEATURES_VALIDATION_ERROR=(1<<2),
  BASIC_VALIDATION_ERROR=(1<<3),
  IS2D_VALIDATION_ERROR=(1<<4),
  LAYOUT2D_VALIDATION_ERROR=(1<<5),
  STEREO_VALIDATION_ERROR=(1<<6),
  VALIDATION_ERROR=(
    FEATURES_VALIDATION_ERROR
    | BASIC_VALIDATION_ERROR
    | IS2D_VALIDATION_ERROR
    | LAYOUT2D_VALIDATION_ERROR
    | STEREO_VALIDATION_ERROR
  ),
  PREPARE_FOR_STANDARDIZATION_ERROR=(1<<7),
  METAL_STANDARDIZATION_ERROR=(1<<8),
  NORMALIZER_STANDARDIZATION_ERROR=(1<<9),
  FRAGMENT_STANDARDIZATION_ERROR=(1<<10),
  CHARGE_STANDARDIZATION_ERROR=(1<<11),
  STANDARDIZATION_ERROR=(
    METAL_STANDARDIZATION_ERROR
    | NORMALIZER_STANDARDIZATION_ERROR
    | FRAGMENT_STANDARDIZATION_ERROR
    | CHARGE_STANDARDIZATION_ERROR
  ),
  OUTPUT_ERROR=(1<<12),
  PIPELINE_ERROR=(
    INPUT_ERROR
    | PREPARE_FOR_VALIDATION_ERROR
    | VALIDATION_ERROR
    | PREPARE_FOR_STANDARDIZATION_ERROR
    | STANDARDIZATION_ERROR
    | OUTPUT_ERROR
  ),
  SANITIZATION_APPLIED=(1<<23),
  METALS_DISCONNECTED=(1<<24),
  NORMALIZATION_APPLIED=(1<<25),
  FRAGMENTS_REMOVED=(1<<26),
  PROTONATION_CHANGED=(1<<27),
  STRUCTURE_MODIFICATION=(
    SANITIZATION_APPLIED
    | METALS_DISCONNECTED
    | NORMALIZATION_APPLIED
    | FRAGMENTS_REMOVED
    | PROTONATION_CHANGED
  )
};

enum RDKIT_MOLSTANDARDIZE_EXPORT PipelineStage {
  PARSING_INPUT=0,
  PREPARE_FOR_VALIDATION=1,
  VALIDATION=2,
  PREPARE_FOR_STANDARDIZATION=3,
  STANDARDIZATION=4,
  SERIALIZING_OUTPUT=5,
  COMPLETED=6
};

struct RDKIT_MOLSTANDARDIZE_EXPORT PipelineLogEntry {
  PipelineStatus status;
  std::string detail;
};

using PipelineLog = std::vector<PipelineLogEntry>;

struct RDKIT_MOLSTANDARDIZE_EXPORT PipelineResult {
  PipelineStatus status;
  PipelineStage stage;
  PipelineLog log;
  std::string inputMolBlock;
  std::string outputMolBlock;
  std::string parentMolBlock;

  void append(PipelineStatus newStatus, const std::string & info); 
};

class RDKIT_MOLSTANDARDIZE_EXPORT Pipeline {
 private:
  PipelineOptions options;

 public:
  Pipeline() = default;
  explicit Pipeline(const PipelineOptions & o) : options(o) {};

  PipelineResult run(const std::string & molblock) const;

 private:
  RWMOL_SPTR parse(const std::string & molblock, PipelineResult & result) const;
  RWMOL_SPTR prepareForValidation(RWMOL_SPTR mol, PipelineResult & result) const;
  RWMOL_SPTR validate(RWMOL_SPTR mol, PipelineResult & result) const;
  RWMOL_SPTR prepareForStandardization(RWMOL_SPTR mol, PipelineResult & result) const;
  RWMOL_SPTR standardize(RWMOL_SPTR mol, PipelineResult & result) const;
  using RWMOL_SPTR_PAIR = std::pair<RWMOL_SPTR, RWMOL_SPTR>;
  RWMOL_SPTR_PAIR makeParent(RWMOL_SPTR mol, PipelineResult & result) const;
  void serialize(RWMOL_SPTR_PAIR output, PipelineResult & result) const;
};

}  // namespace MolStandardize
}  // namespace RDKit

#endif
