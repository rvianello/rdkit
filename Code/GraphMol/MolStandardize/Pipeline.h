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
  bool allowEnhancedStereo {false};
  double is2DZeroThreshold {1e-3};
  double atomClashLimit {0.3};
  double bondLengthLimit {10.};

  // cleanup/standardization
  std::string metalNof {"[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra]~[#7,#8,F]"};
  std::string metalNon {};

  // serialization
  bool ouputV2000 {false};
};

enum RDKIT_MOLSTANDARDIZE_EXPORT PipelineStatus {
  NO_ERROR=0,
  INPUT_ERROR=(1<<0),
  FEATURES_VALIDATION_ERROR=(1<<1),
  BASIC_VALIDATION_ERROR=(1<<2),
  IS2D_VALIDATION_ERROR=(1<<3),
  LAYOUT2D_VALIDATION_ERROR=(1<<4),
  STEREO_VALIDATION_ERROR=(1<<5),
  VALIDATION_ERROR=(
    FEATURES_VALIDATION_ERROR
    | BASIC_VALIDATION_ERROR
    | IS2D_VALIDATION_ERROR
    | LAYOUT2D_VALIDATION_ERROR
    | STEREO_VALIDATION_ERROR
  ),
  PREPARE_STANDARDIZATION_ERROR=(1<<6),
  METAL_STANDARDIZATION_ERROR=(1<<7),
  NORMALIZER_STANDARDIZATION_ERROR=(1<<8),
  FRAGMENT_STANDARDIZATION_ERROR=(1<<9),
  CHARGE_STANDARDIZATION_ERROR=(1<<10),
  PROTONATION_STANDARDIZATION_ERROR=(1<<11),
  STANDARDIZATION_ERROR=(
    PREPARE_STANDARDIZATION_ERROR
    | METAL_STANDARDIZATION_ERROR
    | NORMALIZER_STANDARDIZATION_ERROR
    | FRAGMENT_STANDARDIZATION_ERROR
    | CHARGE_STANDARDIZATION_ERROR
    | PROTONATION_STANDARDIZATION_ERROR
  ),
  OUTPUT_ERROR=(1<<12)
};

enum RDKIT_MOLSTANDARDIZE_EXPORT PipelineStage {
  STARTED=0,
  PARSING_INPUT=1,
  VALIDATION=2,
  STANDARDIZATION=3,
  SERIALIZING_OUTPUT=4,
  COMPLETED=5
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

  PipelineResult run(const std::string & molblock);

 private:
  RWMOL_SPTR parse(const std::string & molblock, PipelineResult & result);
  void validate(const ROMol & mol, PipelineResult & result);
  using RWMOL_SPTR_PAIR = std::pair<RWMOL_SPTR, RWMOL_SPTR>;
  RWMOL_SPTR_PAIR standardize(RWMOL_SPTR mol, PipelineResult & result);
  void serialize(RWMOL_SPTR_PAIR output, PipelineResult & result);
};

}  // namespace MolStandardize
}  // namespace RDKit

#endif
