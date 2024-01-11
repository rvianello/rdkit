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
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace RDKit {

class ROMol;
class RWMol;

namespace MolStandardize {

struct RDKIT_MOLSTANDARDIZE_EXPORT PipelineOptions {
  // parsing
  bool strictParsing {false};

  // validation
  bool reportAllFailures {true};
  double is2DZeroThreshold {1e-3};
  double atomClashLimit {0.3};
  double bondLengthLimit {10.};

  // cleanup/standardization

  // serialization
  bool ouputV2000 {false};
};

enum RDKIT_MOLSTANDARDIZE_EXPORT PipelineStatus {
  NO_ERROR=0,
  INPUT_ERROR=(1<<0),
  BASIC_VALIDATION_ERROR=(1<<1),
  IS2D_VALIDATION_ERROR=(1<<2),
  LAYOUT2D_VALIDATION_ERROR=(1<<3),
  STEREO_VALIDATION_ERROR=(1<<4),
  VALIDATION_ERROR=(
    BASIC_VALIDATION_ERROR
    | IS2D_VALIDATION_ERROR
    | LAYOUT2D_VALIDATION_ERROR
    | STEREO_VALIDATION_ERROR
  ),
  PREPARE_STANDARDIZATION_ERROR=(1<<5),
  METAL_STANDARDIZATION_ERROR=(1<<6),
  NORMALIZER_STANDARDIZATION_ERROR=(1<<7),
  PROTONATION_STANDARDIZATION_ERROR=(1<<8),
  FRAGMENT_STANDARDIZATION_ERROR=(1<<9),
  CHARGE_STANDARDIZATION_ERROR=(1<<10),
  STANDARDIZATION_ERROR=(
    PREPARE_STANDARDIZATION_ERROR
    | METAL_STANDARDIZATION_ERROR
    | NORMALIZER_STANDARDIZATION_ERROR
    | PROTONATION_STANDARDIZATION_ERROR
    | FRAGMENT_STANDARDIZATION_ERROR
    | CHARGE_STANDARDIZATION_ERROR
  ),
  OUTPUT_ERROR=(1<<11)
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
  std::unique_ptr<RWMol> parse(const std::string & molblock, PipelineResult & result);
  void validate(const ROMol & mol, PipelineResult & result);
  std::unique_ptr<RWMol> standardize(RWMol & mol, PipelineResult & result);
  void serialize(const ROMol & mol, PipelineResult & result);
};

}  // namespace MolStandardize
}  // namespace RDKit

#endif
