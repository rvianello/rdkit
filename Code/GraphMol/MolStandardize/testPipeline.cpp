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
#include <GraphMol/Chirality.h>
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

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::PARSING_INPUT);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
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

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE((result.status & MolStandardize::STANDARDIZATION_ERROR) == MolStandardize::NORMALIZER_STANDARDIZATION_ERROR);
    REQUIRE(result.status == (
      MolStandardize::BASIC_VALIDATION_ERROR |
      MolStandardize::PREPARE_FOR_STANDARDIZATION_ERROR |
      MolStandardize::NORMALIZER_STANDARDIZATION_ERROR
      ));
  }

  SECTION("failing Isotopes validation") {
    const char * molblock = R"(
  Mrv2311 01312409582D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -15.3955 7.6033 0 0 MASS=3
M  V30 2 C -16.7292 6.8333 0 0
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::BASIC_VALIDATION_ERROR);
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::FEATURES_VALIDATION_ERROR);
  }

  SECTION("failing features validation, aromatic bonds") {
    const char * molblock = R"(
  Mrv2311 02272411562D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -10.3542 4.29 0 0
M  V30 2 C -11.6879 3.52 0 0
M  V30 3 C -11.6879 1.9798 0 0
M  V30 4 N -10.3542 1.21 0 0
M  V30 5 C -9.0204 1.9798 0 0
M  V30 6 C -9.0204 3.52 0 0
M  V30 7 C -10.3542 5.83 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 1 6
M  V30 3 4 2 3
M  V30 4 4 5 6
M  V30 5 1 1 7
M  V30 6 4 3 4
M  V30 7 4 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::FEATURES_VALIDATION_ERROR);
  }

  SECTION("failing radical validation, disallowed radical") {
    const char * molblock = R"(
  Mrv2311 02082417212D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -20.9372 7.145 0 0 RAD=2
M  V30 2 C -22.2708 6.375 0 0
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::BASIC_VALIDATION_ERROR);
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
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
M  V30 1 C -3.05 5.48 0 0
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::LAYOUT2D_VALIDATION_ERROR);
  }

  SECTION("failing validation, abnormally long bonds") {
    const char * molblock = R"(
          01112413352D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -28.1663 10.4367 0 0
M  V30 2 C -29.5 9.6667 0 0
M  V30 3 C -29.5 11.2067 0 0
M  V30 4 F 150.0 10.4367 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3
M  V30 3 1 3 1
M  V30 4 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::LAYOUT2D_VALIDATION_ERROR);

    molblock = R"(
  Mrv2311 02222409302D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 17 17 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -9.0205 2.1033 0 0
M  V30 2 C -10.3542 1.3333 0 0
M  V30 3 C -7.4805 2.1033 0 0
M  V30 4 C -5.9405 2.1033 0 0
M  V30 5 C -4.4005 2.1033 0 0
M  V30 6 C -2.8605 2.1033 0 0
M  V30 7 C -1.3205 2.1033 0 0
M  V30 8 C 0.2195 2.1033 0 0
M  V30 9 C 1.7595 2.1033 0 0
M  V30 10 C 3.2995 2.1033 0 0
M  V30 11 C 4.8395 2.1033 0 0
M  V30 12 C 6.3795 2.1033 0 0
M  V30 13 C 7.9195 2.1033 0 0
M  V30 14 C 9.4595 2.1033 0 0
M  V30 15 C 10.9995 2.1033 0 0
M  V30 16 C 12.5395 2.1033 0 0
M  V30 17 C 13.7854 1.1981 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 13 14
M  V30 14 1 14 15
M  V30 15 1 15 16
M  V30 16 1 16 17
M  V30 17 1 17 2
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    // long bonds in rings are by defaul allowed
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
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
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::IS2D_VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::STEREO_VALIDATION_ERROR);
  }
}

TEST_CASE("VALIDATION_WITH_ALLOW_EMPTY_MOLS_OPTION") {

  MolStandardize::PipelineOptions options;
  options.allowEmptyMolecules = true;
  MolStandardize::Pipeline pipeline(options);

  SECTION("no atoms produces no error") {
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
    REQUIRE(result.status == MolStandardize::NO_EVENT);
  }
}

TEST_CASE("VALIDATION_WITH_DISALLOWED_LONG_BONDS_IN_RINGS") {

  MolStandardize::PipelineOptions options;
  options.bondLengthLimit = 10.; // adapted to test structure
  options.allowLongBondsInRings = false;
  MolStandardize::Pipeline pipeline(options);

  SECTION("report long bonds in rings") {
    const char * molblock = R"(
  Mrv2311 02222409302D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 17 17 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -9.0205 2.1033 0 0
M  V30 2 C -10.3542 1.3333 0 0
M  V30 3 C -7.4805 2.1033 0 0
M  V30 4 C -5.9405 2.1033 0 0
M  V30 5 C -4.4005 2.1033 0 0
M  V30 6 C -2.8605 2.1033 0 0
M  V30 7 C -1.3205 2.1033 0 0
M  V30 8 C 0.2195 2.1033 0 0
M  V30 9 C 1.7595 2.1033 0 0
M  V30 10 C 3.2995 2.1033 0 0
M  V30 11 C 4.8395 2.1033 0 0
M  V30 12 C 6.3795 2.1033 0 0
M  V30 13 C 7.9195 2.1033 0 0
M  V30 14 C 9.4595 2.1033 0 0
M  V30 15 C 10.9995 2.1033 0 0
M  V30 16 C 12.5395 2.1033 0 0
M  V30 17 C 13.7854 1.1981 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 13 14
M  V30 14 1 14 15
M  V30 15 1 15 16
M  V30 16 1 16 17
M  V30 17 1 17 2
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::VALIDATION_ERROR);
    REQUIRE(result.status & MolStandardize::LAYOUT2D_VALIDATION_ERROR);
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

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::METALS_DISCONNECTED);

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

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    // nitro groups are sanitized in a pre-validation step.
    // this test case is not expected to register any errors.
    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::SANITIZATION_APPLIED);

    std::unique_ptr<RWMol> mol(MolBlockToMol(result.outputMolBlock, false, false));
    REQUIRE(mol);

    std::string smiles {MolToSmiles(*mol)};
    REQUIRE(smiles == "C[N+](=O)[O-]");
  }

  SECTION("Phosphate normalization") {

    const char * molblock_a = R"(
  Mrv2311 04152413292D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 P -15.6247 3.9575 0 0 CHG=1
M  V30 2 C -16.9583 3.1875 0 0
M  V30 3 O -15.6247 5.4975 0 0 CHG=-1
M  V30 4 S -14.291 3.1875 0 0 CHG=-1
M  V30 5 C -15.6247 2.4175 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 5
M  V30 3 1 1 3
M  V30 4 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result_a = pipeline.run(molblock_a);

    for (auto & info : result_a.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    // nitro groups are sanitized in a pre-validation step.
    // this test case is not expected to register any errors.
    REQUIRE(result_a.stage == MolStandardize::COMPLETED);
    REQUIRE((result_a.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result_a.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(result_a.status & MolStandardize::NORMALIZATION_APPLIED);
    REQUIRE(result_a.status & MolStandardize::PROTONATION_CHANGED);

    std::unique_ptr<RWMol> mol_a(MolBlockToMol(result_a.outputMolBlock, false, false));
    REQUIRE(mol_a);

    std::string smiles_a {MolToSmiles(*mol_a)};
    REQUIRE(smiles_a == "CP(C)(=O)S");

    const char * molblock_b = R"(
  Mrv2311 04152413292D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 P -15.6247 3.9575 0 0 CHG=1
M  V30 2 C -16.9583 3.1875 0 0
M  V30 3 S -15.6247 5.4975 0 0 CHG=-1
M  V30 4 O -14.291 3.1875 0 0 CHG=-1
M  V30 5 C -15.6247 2.4175 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 5
M  V30 3 1 1 3
M  V30 4 1 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result_b = pipeline.run(molblock_b);

    for (auto & info : result_b.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    // nitro groups are sanitized in a pre-validation step.
    // this test case is not expected to register any errors.
    REQUIRE(result_b.stage == MolStandardize::COMPLETED);
    REQUIRE((result_b.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result_b.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(result_b.status & MolStandardize::NORMALIZATION_APPLIED);
    REQUIRE(result_b.status & MolStandardize::PROTONATION_CHANGED);

    std::unique_ptr<RWMol> mol_b(MolBlockToMol(result_b.outputMolBlock, false, false));
    REQUIRE(mol_b);

    std::string smiles_b {MolToSmiles(*mol_b)};
    REQUIRE(smiles_b == "CP(C)(=O)S");
  }

  SECTION("normalize w/ RDKit's default normalizer transformations") {
    const char * molblock = R"(
  Mrv2311 03112410152D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 S -20.083 2.4575 0 0
M  V30 2 C -21.4167 1.6875 0 0
M  V30 3 O -20.083 3.9975 0 0
M  V30 4 C -18.7493 1.6875 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 4
M  V30 3 2 1 3
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineOptions options;
    options.normalizerData = "";
    MolStandardize::Pipeline customizedPipeline(options);

    MolStandardize::PipelineResult result = customizedPipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(result.status & MolStandardize::NORMALIZATION_APPLIED);

    std::unique_ptr<RWMol> mol(MolBlockToMol(result.outputMolBlock, false, false));
    REQUIRE(mol);

    std::string smiles {MolToSmiles(*mol)};
    REQUIRE(smiles == "C[S+](C)[O-]");
  }

  SECTION("normalization of 1,3- 1,5- conjugated systems favors application within rings") {
    const char * molblock {};
    MolStandardize::PipelineResult result;
    std::unique_ptr<RWMol> parentMol;
    std::string parentSmiles;

    // 1,3- conjugated cation - test a first ctab permutation
    molblock = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -11.833300 8.935800 0.000000 0 CHG=1
M  V30 2 C -13.167100 8.165800 0.000000 0
M  V30 3 C -11.833300 5.855800 0.000000 0
M  V30 4 C -10.499600 6.625600 0.000000 0
M  V30 5 C -10.499600 8.165800 0.000000 0
M  V30 6 C -11.833300 10.475800 0.000000 0
M  V30 7 N -13.167100 6.625600 0.000000 0
M  V30 8 N -14.500700 8.935800 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 4
M  V30 2 1 4 5
M  V30 3 2 1 2
M  V30 4 1 1 5
M  V30 5 1 1 6
M  V30 6 1 2 8
M  V30 7 1 2 7
M  V30 8 1 7 3
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (MolStandardize::NORMALIZATION_APPLIED | MolStandardize::PROTONATION_CHANGED)
    );

    parentMol.reset(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    parentSmiles = MolToSmiles(*parentMol);
    REQUIRE(parentSmiles == "CN1CCCN=C1N");

    // 1,3- conjugated cation - test a second ctab permutation
    molblock = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -11.833300 8.935800 0.000000 0 CHG=1
M  V30 2 C -13.167100 8.165800 0.000000 0
M  V30 3 C -11.833300 5.855800 0.000000 0
M  V30 4 C -10.499600 6.625600 0.000000 0
M  V30 5 C -10.499600 8.165800 0.000000 0
M  V30 6 C -11.833300 10.475800 0.000000 0
M  V30 7 N -14.500700 8.935800 0.000000 0
M  V30 8 N -13.167100 6.625600 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 4
M  V30 2 1 4 5
M  V30 3 2 1 2
M  V30 4 1 1 5
M  V30 5 1 1 6
M  V30 6 1 2 7
M  V30 7 1 2 8
M  V30 8 1 8 3
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (MolStandardize::NORMALIZATION_APPLIED | MolStandardize::PROTONATION_CHANGED)
    );

    parentMol.reset(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    parentSmiles = MolToSmiles(*parentMol);
    REQUIRE(parentSmiles == "CN1CCCN=C1N");

    // 1,5- conjugated cation - test a first ctab permutation
    molblock = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -3.958300 0.790000 0.000000 0 CHG=1
M  V30 2 C -5.292100 0.020000 0.000000 0
M  V30 3 C -5.292100 -1.520200 0.000000 0
M  V30 4 C -3.958300 -2.290000 0.000000 0
M  V30 5 C -2.624600 0.020000 0.000000 0
M  V30 6 N -2.624600 -1.520200 0.000000 0
M  V30 7 C -3.958300 2.330000 0.000000 0
M  V30 8 N -3.958300 -3.830000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 2 3 4
M  V30 3 2 1 2
M  V30 4 1 1 5
M  V30 5 1 1 7
M  V30 6 1 4 8
M  V30 7 1 6 5
M  V30 8 1 4 6
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (MolStandardize::NORMALIZATION_APPLIED | MolStandardize::PROTONATION_CHANGED)
    );

    parentMol.reset(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    parentSmiles = MolToSmiles(*parentMol);
    REQUIRE(parentSmiles == "CN1C=CC(N)=NC1");

    // 1,5- conjugated cation - test a second ctab permutation
    molblock = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -3.958300 0.790000 0.000000 0 CHG=1
M  V30 2 C -5.292100 0.020000 0.000000 0
M  V30 3 C -5.292100 -1.520200 0.000000 0
M  V30 4 C -3.958300 -2.290000 0.000000 0
M  V30 5 C -2.624600 0.020000 0.000000 0
M  V30 6 N -3.958300 -3.830000 0.000000 0
M  V30 7 N -2.624600 -1.520200 0.000000 0
M  V30 8 C -3.958300 2.330000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 2 3 4
M  V30 3 2 1 2
M  V30 4 1 1 5
M  V30 5 1 1 8
M  V30 6 1 4 6
M  V30 7 1 7 5
M  V30 8 1 4 7
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (MolStandardize::NORMALIZATION_APPLIED | MolStandardize::PROTONATION_CHANGED)
    );

    parentMol.reset(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    parentSmiles = MolToSmiles(*parentMol);
    REQUIRE(parentSmiles == "CN1C=CC(N)=NC1");
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

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (MolStandardize::METALS_DISCONNECTED | MolStandardize::FRAGMENTS_REMOVED)
    );

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "NCC(=O)O");

    std::unique_ptr<RWMol> outputMol(MolBlockToMol(result.outputMolBlock, false, false));
    REQUIRE(outputMol);
    std::string outputSmiles {MolToSmiles(*outputMol)};
    REQUIRE(outputSmiles == "[NH3+]CC(=O)[O-]");
  }

  SECTION("standardize zwitterion with quaternary nitrogen") {
    const char * molblock = R"(
  Mrv2311 02052411472D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -10.9997 4.77 0 0 CHG=1
M  V30 2 C -12.3333 4 0 0
M  V30 3 C -10.9997 6.31 0 0
M  V30 4 C -9.666 4 0 0
M  V30 5 C -10.9997 3.23 0 0
M  V30 6 C -9.666 2.46 0 0
M  V30 7 O -10.9997 1.69 0 0
M  V30 8 O -8.3323 1.69 0 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 4 6
M  V30 2 2 6 7
M  V30 3 1 6 8
M  V30 4 1 2 1
M  V30 5 1 1 3
M  V30 6 1 1 4
M  V30 7 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) == MolStandardize::NO_EVENT);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "C[N+](C)(C)CC(=O)O");

    std::unique_ptr<RWMol> outputMol(MolBlockToMol(result.outputMolBlock, false, false));
    REQUIRE(outputMol);
    std::string outputSmiles {MolToSmiles(*outputMol)};
    REQUIRE(outputSmiles == "C[N+](C)(C)CC(=O)[O-]");
  }

  SECTION("uncharge tertiary amine w/ explicit hydrogen") {
    const char * molblock = R"(
  Mrv2311 02012412352D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -13.2705 5.77 0 0
M  V30 2 N -14.6042 5 0 0 CHG=1
M  V30 3 H -15.9378 5.77 0 0
M  V30 4 C -14.6042 3.46 0 0
M  V30 5 C -13.2705 4.23 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 4
M  V30 3 1 2 5
M  V30 4 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) == MolStandardize::PROTONATION_CHANGED);
    REQUIRE(result.parentMolBlock == result.outputMolBlock);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "CN(C)C");
  }

  SECTION("standardize preserves explicit Hs on chiral centers") {
    const char * molblock = R"(
  Mrv2311 03112410582D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 19 19 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -22.6873 4.415 0 0
M  V30 2 C -24.0208 3.6451 0 0
M  V30 3 C -24.0208 2.1049 0 0
M  V30 4 C -22.6873 1.3349 0 0
M  V30 5 C -21.3536 2.1049 0 0
M  V30 6 C -21.3536 3.6451 0 0
M  V30 7 S -20.02 4.415 0 0 CHG=1
M  V30 8 O -20.02 5.955 0 0 CHG=-1
M  V30 9 C -18.6863 3.6451 0 0
M  V30 10 C -17.3526 4.415 0 0
M  V30 11 C -16.0189 3.6451 0 0
M  V30 12 C -14.6853 4.4151 0 0 CFG=1
M  V30 13 C -13.3516 3.6452 0 0
M  V30 14 F -14.6853 5.9551 0 0
M  V30 15 H -14.6853 2.8752 0 0
M  V30 16 C -25.3545 1.3349 0 0
M  V30 17 O -26.6882 2.105 0 0
M  V30 18 O -25.3546 -0.2051 0 0
M  V30 19 Na -28.0219 1.335 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 6 7
M  V30 8 1 7 9
M  V30 9 1 7 8
M  V30 10 1 9 10
M  V30 11 1 10 11
M  V30 12 1 11 12
M  V30 13 1 12 13
M  V30 14 1 12 15 CFG=1
M  V30 15 1 12 14
M  V30 16 1 3 16
M  V30 17 1 16 17
M  V30 18 2 16 18
M  V30 19 1 17 19
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (
        MolStandardize::METALS_DISCONNECTED
        | MolStandardize::FRAGMENTS_REMOVED
        | MolStandardize::PROTONATION_CHANGED
        )
    );
    REQUIRE(result.outputMolBlock == result.parentMolBlock);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "[H][C@](C)(F)CCC[S+]([O-])C1=CC=C(C(=O)O)C=C1");
  }

  SECTION("standardize preserves isotopically marked explicit Hs") {
    const char * molblock = R"(
  Mrv2311 03112410572D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 17 17 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -22.6873 4.415 0 0
M  V30 2 C -24.0208 3.6451 0 0
M  V30 3 C -24.0208 2.1049 0 0
M  V30 4 C -22.6873 1.3349 0 0
M  V30 5 C -21.3536 2.1049 0 0
M  V30 6 C -21.3536 3.6451 0 0
M  V30 7 S -20.02 4.415 0 0 CHG=1
M  V30 8 O -20.02 5.955 0 0 CHG=-1
M  V30 9 C -18.6863 3.6451 0 0
M  V30 10 C -17.3526 4.415 0 0
M  V30 11 C -16.0189 3.6451 0 0
M  V30 12 C -14.6853 4.4151 0 0
M  V30 13 H -13.3516 3.6452 0 0 MASS=2
M  V30 14 C -25.3545 1.3349 0 0
M  V30 15 O -26.6882 2.105 0 0
M  V30 16 O -25.3546 -0.2051 0 0
M  V30 17 Na -28.0219 1.335 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 6 7
M  V30 8 1 7 9
M  V30 9 1 7 8
M  V30 10 1 9 10
M  V30 11 1 10 11
M  V30 12 1 11 12
M  V30 13 1 3 14
M  V30 14 1 14 15
M  V30 15 2 14 16
M  V30 16 1 15 17
M  V30 17 1 12 13
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (
        MolStandardize::METALS_DISCONNECTED
        //| MolStandardize::NORMALIZATION_APPLIED
        | MolStandardize::FRAGMENTS_REMOVED
        | MolStandardize::PROTONATION_CHANGED
        )
    );
    REQUIRE(result.outputMolBlock == result.parentMolBlock);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "[2H]CCCC[S+]([O-])C1=CC=C(C(=O)O)C=C1");
  }

  SECTION("standardize preserves generic explicit Hs") {
    const char * molblock = R"(
  Mrv2311 03112410542D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 17 17 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -22.6873 4.415 0 0
M  V30 2 C -24.0208 3.6451 0 0
M  V30 3 C -24.0208 2.1049 0 0
M  V30 4 C -22.6873 1.3349 0 0
M  V30 5 C -21.3536 2.1049 0 0
M  V30 6 C -21.3536 3.6451 0 0
M  V30 7 S -20.02 4.415 0 0 CHG=1
M  V30 8 O -20.02 5.955 0 0 CHG=-1
M  V30 9 C -18.6863 3.6451 0 0
M  V30 10 C -17.3526 4.415 0 0
M  V30 11 C -16.0189 3.6451 0 0
M  V30 12 C -14.6853 4.4151 0 0
M  V30 13 H -13.3516 3.6452 0 0
M  V30 14 C -25.3545 1.3349 0 0
M  V30 15 O -26.6882 2.105 0 0
M  V30 16 O -25.3546 -0.2051 0 0
M  V30 17 Na -28.0219 1.335 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 6 7
M  V30 8 1 7 9
M  V30 9 1 7 8
M  V30 10 1 9 10
M  V30 11 1 10 11
M  V30 12 1 11 12
M  V30 13 1 3 14
M  V30 14 1 14 15
M  V30 15 2 14 16
M  V30 16 1 15 17
M  V30 17 1 12 13
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (
        MolStandardize::METALS_DISCONNECTED
        | MolStandardize::FRAGMENTS_REMOVED
        | MolStandardize::PROTONATION_CHANGED
        )
    );
    REQUIRE(result.outputMolBlock == result.parentMolBlock);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "[H]CCCC[S+]([O-])C1=CC=C(C(=O)O)C=C1");
  }

  SECTION("standardize doesn't remove wedged bonds from non-stereogenic centers") {
    const char * molblock = R"(
  Mrv2311 03112410512D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 16 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -22.6875 4.415 0 0
M  V30 2 C -24.021 3.6452 0 0
M  V30 3 C -24.021 2.1049 0 0
M  V30 4 C -22.6875 1.3349 0 0
M  V30 5 C -21.3537 2.1049 0 0
M  V30 6 C -21.3537 3.6452 0 0
M  V30 7 S -20.0202 4.415 0 0 CHG=1
M  V30 8 O -20.0202 5.955 0 0 CHG=-1
M  V30 9 C -18.6865 3.6452 0 0
M  V30 10 C -17.3528 4.415 0 0 CFG=2
M  V30 11 C -16.019 3.6452 0 0
M  V30 12 C -25.3547 1.3349 0 0
M  V30 13 O -26.6885 2.105 0 0
M  V30 14 O -25.3548 -0.2051 0 0
M  V30 15 Na -28.0222 1.3351 0 0
M  V30 16 C -17.3528 5.955 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 6 7
M  V30 8 1 7 9
M  V30 9 1 7 8
M  V30 10 1 9 10
M  V30 11 1 10 11 CFG=1
M  V30 12 1 3 12
M  V30 13 2 12 14
M  V30 14 1 12 13
M  V30 15 1 13 15
M  V30 16 1 10 16
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (
        MolStandardize::METALS_DISCONNECTED
        | MolStandardize::FRAGMENTS_REMOVED
        | MolStandardize::PROTONATION_CHANGED
        )
    );
    REQUIRE(result.outputMolBlock == result.parentMolBlock);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "CC(C)C[S+]([O-])C1=CC=C(C(=O)O)C=C1");
  
    Chirality::reapplyMolBlockWedging(*parentMol);

    const Bond * wedged = nullptr;
    for (auto bond: parentMol->bonds()) {
      auto bondDir = bond->getBondDir();
      if (bondDir == Bond::BondDir::BEGINWEDGE) {
        wedged = bond;
        break;
      }
    }
    REQUIRE(wedged != nullptr);

    auto beginAtom = wedged->getBeginAtom();
    REQUIRE(beginAtom->getAtomicNum() == 6);
    REQUIRE(beginAtom->getDegree() == 3);
    auto endAtom = wedged->getEndAtom();
    REQUIRE(endAtom->getAtomicNum() == 6);
    REQUIRE(endAtom->getDegree() == 1);
  }

  SECTION("standardize doesn't remove wavy bonds from non-stereogenic centers") {
    const char * molblock = R"(
  Mrv2311 03112410452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 16 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -22.6875 4.415 0 0
M  V30 2 C -24.021 3.6452 0 0
M  V30 3 C -24.021 2.1049 0 0
M  V30 4 C -22.6875 1.3349 0 0
M  V30 5 C -21.3537 2.1049 0 0
M  V30 6 C -21.3537 3.6452 0 0
M  V30 7 S -20.0202 4.415 0 0 CHG=1
M  V30 8 O -20.0202 5.955 0 0 CHG=-1
M  V30 9 C -18.6865 3.6452 0 0
M  V30 10 C -17.3528 4.415 0 0 CFG=3
M  V30 11 C -16.019 3.6452 0 0
M  V30 12 C -25.3547 1.3349 0 0
M  V30 13 O -26.6885 2.105 0 0
M  V30 14 O -25.3548 -0.2051 0 0
M  V30 15 Na -28.0222 1.3351 0 0
M  V30 16 C -17.3528 5.955 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 6 7
M  V30 8 1 7 9
M  V30 9 1 7 8
M  V30 10 1 9 10
M  V30 11 1 10 11 CFG=2
M  V30 12 1 3 12
M  V30 13 2 12 14
M  V30 14 1 12 13
M  V30 15 1 13 15
M  V30 16 1 10 16
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION)
      == (
        MolStandardize::METALS_DISCONNECTED
        | MolStandardize::FRAGMENTS_REMOVED
        | MolStandardize::PROTONATION_CHANGED
        )
    );
    REQUIRE(result.outputMolBlock == result.parentMolBlock);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "CC(C)C[S+]([O-])C1=CC=C(C(=O)O)C=C1");
  
    Chirality::reapplyMolBlockWedging(*parentMol);

    const Bond * wedged = nullptr;
    for (auto bond: parentMol->bonds()) {
      auto bondDir = bond->getBondDir();
      if (bondDir == Bond::BondDir::UNKNOWN) {
        wedged = bond;
        break;
      }
    }
    REQUIRE(wedged != nullptr);

    auto beginAtom = wedged->getBeginAtom();
    REQUIRE(beginAtom->getAtomicNum() == 6);
    REQUIRE(beginAtom->getDegree() == 3);
    auto endAtom = wedged->getEndAtom();
    REQUIRE(endAtom->getAtomicNum() == 6);
    REQUIRE(endAtom->getDegree() == 1);
  }

  SECTION("standardize removes wavy bonds from double bonds w/ stereo type 'either'") {
    const char * molblock = R"(
  Mrv2311 04172413232D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -13.9163 3.9158 0 0
M  V30 2 C -15.25 3.1458 0 0
M  V30 3 C -13.9163 5.4558 0 0
M  V30 4 C -15.25 6.2258 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 2 1 3
M  V30 3 1 3 4 CFG=2
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION) == MolStandardize::NORMALIZATION_APPLIED
    );
    REQUIRE(result.outputMolBlock == result.parentMolBlock);

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "CC=CC");
  
    Chirality::reapplyMolBlockWedging(*parentMol);

    // no wavy bond is expected to be found
    const Bond * wavy = nullptr;
    for (auto bond: parentMol->bonds()) {
      auto bondDir = bond->getBondDir();
      if (bondDir == Bond::BondDir::UNKNOWN) {
        wavy = bond;
        break;
      }
    }
    REQUIRE(wavy == nullptr);

    // the double bond should have stereo type STEREOANY
    const Bond * doubleBond = nullptr;
    for (auto bond: parentMol->bonds()) {
      auto bondType = bond->getBondType();
      if (bondType == Bond::DOUBLE) {
        doubleBond = bond;
        break;
      }
    }
    REQUIRE(doubleBond != nullptr);
    REQUIRE(doubleBond->getStereo() == Bond::STEREOANY);
  }

  SECTION("pipeline doesn't remove stereo bonds from biaryls") {
    const char * molblock = R"(
  Mrv2311 02092409022D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 13 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -17.8543 8.7068 0 0
M  V30 2 C -19.1878 7.9368 0 0
M  V30 3 C -19.1878 6.3966 0 0
M  V30 4 C -17.8543 5.6266 0 0
M  V30 5 C -16.5205 6.3966 0 0
M  V30 6 C -16.5205 7.9368 0 0
M  V30 7 C -17.8543 4.0866 0 0
M  V30 8 C -19.1879 3.3166 0 0
M  V30 9 C -19.1879 1.7764 0 0
M  V30 10 C -17.8544 1.0064 0 0
M  V30 11 C -16.5206 1.7763 0 0
M  V30 12 C -16.5206 3.3165 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5 CFG=1
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 4 7
M  V30 8 1 8 9
M  V30 9 2 9 10
M  V30 10 1 10 11
M  V30 11 2 11 12
M  V30 12 2 7 8
M  V30 13 1 12 7
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) == MolStandardize::NO_EVENT);
    REQUIRE(result.outputMolBlock == result.parentMolBlock);
    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
  
    Chirality::reapplyMolBlockWedging(*parentMol);

    const Bond * wedged = nullptr;
    for (auto bond: parentMol->bonds()) {
      auto bondDir = bond->getBondDir();
      if (bondDir == Bond::BondDir::BEGINWEDGE) {
        wedged = bond;
        break;
      }
    }
    REQUIRE(wedged != nullptr);

    auto beginAtom = wedged->getBeginAtom();
    // there's only two position with degree 3
    // and they are equivalent
    REQUIRE(beginAtom->getDegree() == 3);
  }

  SECTION("SO2H normalization") {
    const char * molblock = R"(
  Mrv2311 03122408072D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 S -10.708 4.2075 0 0
M  V30 2 C -12.0417 3.4375 0 0
M  V30 3 O -10.708 5.7475 0 0
M  V30 4 O -9.3743 3.4375 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 2 1 3
M  V30 3 2 1 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) != MolStandardize::NO_EVENT);
    REQUIRE(
      (result.status & MolStandardize::STRUCTURE_MODIFICATION) == MolStandardize::NORMALIZATION_APPLIED);

    REQUIRE(result.outputMolBlock == result.parentMolBlock);

    std::unique_ptr<RWMol> inputMol(MolBlockToMol(result.inputMolBlock, false, false));
    REQUIRE(inputMol);
    std::string inputSmiles {MolToSmiles(*inputMol)};
    REQUIRE(inputSmiles == "C[SH](=O)=O");

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "CS(=O)O");
  }

  SECTION("Partial disconnection of metals") {

    // This test requires handling the disconnection of
    // bivalent metals
    MolStandardize::PipelineOptions options;
    options.metalNof = "[Li,Na,K,Be,Mg,Ca]~[#7,#8,F]";
    MolStandardize::Pipeline pipeline(options);


    const char * molblock = R"(
  Mrv2311 03182415572D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 16 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.3342 -4.235 0 0
M  V30 2 C -3.3342 -2.695 0 0
M  V30 3 O -2.0005 -1.925 0 0
M  V30 4 Mg -2.0005 -0.385 0 0
M  V30 5 C -2.0005 1.155 0 0
M  V30 6 C -3.3342 1.925 0 0
M  V30 7 C -3.3342 3.465 0 0
M  V30 8 C -2.0005 4.235 0 0
M  V30 9 C -0.6668 3.465 0 0
M  V30 10 C -0.6668 1.925 0 0
M  V30 11 C 0.6668 1.155 0 0
M  V30 12 C 0.6668 -0.385 0 0
M  V30 13 C 2.0005 -1.155 0 0
M  V30 14 O 2.0005 -2.695 0 0
M  V30 15 O 3.3342 -0.385 0 0
M  V30 16 Na 3.3342 1.155 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 5 10
M  V30 11 1 10 11
M  V30 12 1 11 12
M  V30 13 1 12 13
M  V30 14 2 13 14
M  V30 15 1 13 15
M  V30 16 1 15 16
M  V30 END BOND
M  V30 END CTAB
M  END
)";

    MolStandardize::PipelineResult result = pipeline.run(molblock);

    for (auto & info : result.log) {
      std::cerr << info.status << " " << info.detail << std::endl;
    }

    REQUIRE(result.stage == MolStandardize::COMPLETED);
    REQUIRE((result.status & MolStandardize::PIPELINE_ERROR) == MolStandardize::NO_EVENT);
    REQUIRE((result.status & MolStandardize::STRUCTURE_MODIFICATION) == MolStandardize::NO_EVENT);

    REQUIRE(result.outputMolBlock == result.parentMolBlock);

    std::unique_ptr<RWMol> inputMol(MolBlockToMol(result.inputMolBlock, false, false));
    REQUIRE(inputMol);
    std::string inputSmiles {MolToSmiles(*inputMol)};
    REQUIRE(inputSmiles == "CCO[Mg]C1CCCCC1CCC(=O)O[Na]");

    std::unique_ptr<RWMol> parentMol(MolBlockToMol(result.parentMolBlock, false, false));
    REQUIRE(parentMol);
    std::string parentSmiles {MolToSmiles(*parentMol)};
    REQUIRE(parentSmiles == "CCO[Mg]C1CCCCC1CCC(=O)O[Na]");
  }
}
