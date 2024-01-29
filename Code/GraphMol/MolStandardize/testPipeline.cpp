//
//  Copyright (C) 2023 Novartis Biomedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolStandardize/Pipeline.h>
#include <memory>
#include <string>

using namespace RDKit;

TEST_CASE("VALIDATION_ERRORS_WITH_DEFAULT_OPTIONS") {

  MolStandardize::Pipeline pipeline;

  SECTION("parse error") {
    const char * molblock = R"(
             sldfj;ldskfj sldkjfsd;lkf 
M  V30 BEGIN CTAB
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    REQUIRE(result.stage == MolStandardize::PARSING_INPUT);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::INPUT_ERROR);
  }

  SECTION("failing RDKit validation, no atoms") {
    const char * molblock = R"(
          10052313452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 0 0 0 0 0
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::BASIC_VALIDATION_ERROR);
  }

  SECTION("failing RDKit validation, bad valence status") {
    const char * molblock = R"(
          10242314442D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.6247 7.5825 0 0
M  V30 2 N -2.9583 6.8125 0 0
M  V30 3 C -4.292 7.5825 0 0
M  V30 4 C -2.9583 5.2725 0 0
M  V30 5 C -1.6247 6.0425 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::STANDARDIZATION);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::STANDARDIZATION_ERROR);
    REQUIRE(result.status == (
      MolStandardize::BASIC_VALIDATION_ERROR
      | MolStandardize::FRAGMENT_STANDARDIZATION_ERROR
      ));
  }

  SECTION("failing features validation, query atom") {
    const char * molblock = R"(
  Mrv2311 01162413552D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 R# -17.3747 6.9367 0 0 RGROUPS=(1 0)
M  V30 2 C -18.7083 6.1667 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::FEATURES_VALIDATION_ERROR);
  }

  SECTION("failing features validation, enhanced stereo") {
    const char * molblock = R"(
  Mrv2311 01162411552D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -18.208 8.52 0 0 CFG=2
M  V30 2 F -19.5417 7.75 0 0
M  V30 3 C -16.8743 7.75 0 0
M  V30 4 Cl -18.208 10.06 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 3 CFG=1
M  V30 2 1 2 1
M  V30 3 1 1 4
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STERAC1 ATOMS=(1 1)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::FEATURES_VALIDATION_ERROR);
  }

  SECTION("failing 2D validation, non-null Z coords") {
    const char * molblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8753 4.9367 0 0
M  V30 2 C -0.4583 4.1667 0.2 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::IS2D_VALIDATION_ERROR);
  }

  SECTION("failing validation, clashing atoms") {
    const char * molblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.083 5.4575 0 0
M  V30 2 C -4.4167 4.6875 0 0
M  V30 3 C -4.3289 6.3627 0 0
M  V30 4 C -3.0 5.5 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::LAYOUT2D_VALIDATION_ERROR);
  }

  SECTION("failing stereo validation, too many stereo bonds (3 subst. case)") {
    const char * molblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.583 5.7075 0 0 CFG=1
M  V30 2 C -2.9167 4.9375 0 0
M  V30 3 C -1.583 7.2475 0 0
M  V30 4 C -0.2493 4.9375 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3 CFG=1
M  V30 3 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::STEREO_VALIDATION_ERROR);
  }

  SECTION("failing stereo validation, adjacent stereo bonds w/ like direction (4 subst. case)") {
    const char * molblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.583 5.7075 0 0
M  V30 2 C -2.9167 4.9375 0 0
M  V30 3 C -1.583 7.2475 0 0
M  V30 4 C -0.2493 4.9375 0 0
M  V30 5 C -1.583 4.1675 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3 CFG=1
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::STEREO_VALIDATION_ERROR);
  }

  SECTION("failing validation, not 2D + adjacent stereo bonds w/ like direction (4 subst. case)") {
    const char * molblock = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.583 5.7075 0 0
M  V30 2 C -2.9167 4.9375 0 0
M  V30 3 C -1.583 7.2475 0 0
M  V30 4 C -0.2493 4.9375 0.5 0
M  V30 5 C -1.583 4.1675 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3 CFG=1
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status != MolStandardize::NO_ERROR);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::IS2D_VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::STEREO_VALIDATION_ERROR);
  }
}


TEST_CASE("STANDARDIZATION") {

  MolStandardize::Pipeline pipeline;

  SECTION("disconnect metal") {
    const char * molblock = R"(
          10282320572D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.0413 5.4992 0 0
M  V30 2 C -2.375 4.7292 0 0
M  V30 3 O -1.0413 7.0392 0 0
M  V30 4 O 0.2924 4.7292 0 0
M  V30 5 Na 0.2924 3.1892 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 4
M  V30 3 2 1 3
M  V30 4 1 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::VALIDATION_ERROR) == 0);
    REQUIRE(result.status == MolStandardize::NO_ERROR);

    std::unique_ptr<RWMol> mol(MolBlockToMol(result.outputMolBlock, false, false));
    REQUIRE(mol);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "CC(=O)O");

    std::unique_ptr<RWMol> outputMol(MolBlockToMol(result.outputMolBlock, false, false));
    REQUIRE(outputMol);
    std::string outputSmiles {MolToSmiles(*outputMol)};
    REQUIRE(outputSmiles == "CC(=O)O");
}

  SECTION("normalize nitro") {
    const char * molblock = R"(
          10282320572D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -1.0413 5.4992 0 0
M  V30 2 C -2.375 4.7292 0 0
M  V30 3 O -1.0413 7.0392 0 0
M  V30 4 O 0.2924 4.7292 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 2 1 4
M  V30 3 2 1 3
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status == MolStandardize::BASIC_VALIDATION_ERROR);

    std::unique_ptr<RWMol> mol(MolBlockToMol(result.outputMolBlock, false, false));
    REQUIRE(mol);

    std::string smiles {MolToSmiles(*mol)};
    REQUIRE(smiles == "C[N+](=O)[O-]");
  }

  SECTION("standardize zwitterion") {
    const char * molblock = R"(
          10282320572D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.0413 5.4992 0 0
M  V30 2 C -2.375 4.7292 0 0
M  V30 3 O -1.0413 7.0392 0 0
M  V30 4 O 0.2924 4.7292 0 0
M  V30 5 N -3.7087 5.4992 0 0 CHG=1
M  V30 6 Na 0.2924 3.1892 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 4
M  V30 3 2 1 3
M  V30 4 1 2 5
M  V30 5 1 4 6
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE(result.status == MolStandardize::NO_ERROR);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "NCC(=O)O");

    std::unique_ptr<RWMol> outputMol(MolBlockToMol(result.outputMolBlock, false, false));
    REQUIRE(outputMol);
    std::string outputSmiles {MolToSmiles(*outputMol)};
    REQUIRE(outputSmiles == "[NH3+]CC(=O)[O-]");
  }
}

