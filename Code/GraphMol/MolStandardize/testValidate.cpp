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
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Chirality.h>

#include <iostream>

using namespace RDKit;
using namespace std;
using namespace MolStandardize;

void testRDKitValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing RDKit validation"
                       << std::endl;

  string smi1, smi2, smi3, smi4;
  RDKitValidation vm;

  // testing RDKitDefault
  smi1 = "CO(C)C";
  unique_ptr<ROMol> m1(SmilesToMol(smi1, 0, false));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (auto &query : errout1) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, "
                "is greater than permitted");
  }

  // testing for molecule with no atoms
  smi2 = "";
  unique_ptr<ROMol> m2(SmilesToMol(smi2, 0, false));
  vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
  for (auto &query : errout2) {
    std::string msg = query.what();
    TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  // testing molecule with multiple valency errors
  smi3 = "CO(C)CCN(=O)=O";
  unique_ptr<ROMol> m3(SmilesToMol(smi3, 0, false));
  vector<ValidationErrorInfo> errout3 = vm.validate(*m3, true);
  std::vector<string> msgs1;
  std::vector<string> ans1 = {
      "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is "
      "greater than permitted",
      "INFO: [ValenceValidation] Explicit valence for atom # 5 N, 5, is "
      "greater than permitted"};
  for (auto &query : errout3) {
    msgs1.emplace_back(query.what());
  }
  TEST_ASSERT(msgs1 == ans1);

  // testing molecule with multiple valency errors and only outputting
  // first error
  smi4 = "CO(C)CCN(=O)=O";
  unique_ptr<ROMol> m4(SmilesToMol(smi4, 0, false));
  vector<ValidationErrorInfo> errout4 = vm.validate(*m4, false);
  std::vector<string> msgs2;
  std::vector<string> ans2 = {
      "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is "
      "greater than permitted"};
  for (auto &query : errout4) {
    msgs2.emplace_back(query.what());
  }
  TEST_ASSERT(msgs2 == ans2);
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMolVSValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing MolVS validation"
                       << std::endl;
  string smi1, smi2, smi3, smi4, smi5, smi6;
  MolVSValidation vm;

  // testing MolVSDefault
  // testing for molecule with no atoms
  smi1 = "";
  unique_ptr<ROMol> m1(SmilesToMol(smi1, 0, false));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (auto &query : errout1) {
    std::string msg = query.what();
    TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  smi2 = "O=C([O-])c1ccccc1";
  unique_ptr<ROMol> m2(SmilesToMol(smi2, 0, false));
  vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
  for (auto &query : errout2) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [NeutralValidation] Not an overall neutral system (-1)");
  }

  smi3 = "CN=[NH+]CN=N";
  unique_ptr<ROMol> m3(SmilesToMol(smi3, 0, false));
  vector<ValidationErrorInfo> errout3 = vm.validate(*m3, true);
  for (auto &query : errout3) {
    std::string msg = query.what();
    TEST_ASSERT(
        msg ==
        "INFO: [NeutralValidation] Not an overall neutral system (+1)");  // fix
                                                                          // to
                                                                          // show
                                                                          // +
                                                                          // sign
  }

  smi4 = "[13CH4]";
  unique_ptr<ROMol> m4(SmilesToMol(smi4, 0, false));
  vector<ValidationErrorInfo> errout4 = vm.validate(*m4, true);
  for (auto &query : errout4) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [IsotopeValidation] Molecule contains isotope 13C");
  }

  smi5 = "[2H]C(Cl)(Cl)Cl";
  unique_ptr<ROMol> m5(SmilesToMol(smi5, 0, false));
  vector<ValidationErrorInfo> errout5 = vm.validate(*m5, true);
  for (auto &query : errout5) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [IsotopeValidation] Molecule contains isotope 2H");
  }

  smi6 = "[2H]OC([2H])([2H])[2H]";
  unique_ptr<ROMol> m6(SmilesToMol(smi6, 0, false));
  vector<ValidationErrorInfo> errout6 = vm.validate(*m6, true);
  for (auto &query : errout6) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [IsotopeValidation] Molecule contains isotope 2H");
  }

  std::string smi7 = "COc1cccc(C=N[N-]C(N)=O)c1[O-].O.O.O.O=[U+2]=O";
  unique_ptr<ROMol> m7(SmilesToMol(smi7, 0, false));
  vector<ValidationErrorInfo> errout7 = vm.validate(*m7, true);
  TEST_ASSERT(errout7.size() != 0);
  for (auto &query : errout7) {
    std::string msg = query.what();
    TEST_ASSERT(msg == "INFO: [FragmentValidation] water/hydroxide is present");
  }

  std::string smi8 = "CC(=O)O.NCC(=O)NCCCCCCCCCCNC(=O)CN";
  unique_ptr<ROMol> m8(SmilesToMol(smi8, 0, false));
  vector<ValidationErrorInfo> errout8 = vm.validate(*m8, true);
  TEST_ASSERT(errout8.size() != 0);
  for (auto &query : errout8) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [FragmentValidation] acetate/acetic acid is present");
  }

  std::string smi9 = "N#CC(Br)(Br)C#N.[Br-].[K+]";
  unique_ptr<ROMol> m9(SmilesToMol(smi9, 0, false));
  vector<ValidationErrorInfo> errout9 = vm.validate(*m9, true);
  std::vector<std::string> ans = {
      "INFO: [FragmentValidation] bromine is present",
      "INFO: [FragmentValidation] potassium is present"};
  TEST_ASSERT(errout9.size() == ans.size());
  for (size_t i = 0; i < errout9.size(); ++i) {
    TEST_ASSERT(errout9[i].what() == ans[i]);
  }

  std::string smi10 = "C1COCCO1.O=C(NO)NO";
  unique_ptr<ROMol> m10(SmilesToMol(smi10, 0, false));
  vector<ValidationErrorInfo> errout10 = vm.validate(*m10, true);
  std::vector<std::string> ans10 = {
      "INFO: [FragmentValidation] 1,2-dimethoxyethane is present",
      "INFO: [FragmentValidation] 1,4-dioxane is present"};
  TEST_ASSERT(errout10.size() == ans10.size());
  for (size_t i = 0; i < errout10.size(); ++i) {
    TEST_ASSERT(errout10[i].what() == ans10[i]);
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testMolVSOptions() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing MolVS Options"
                       << std::endl;
  vector<boost::shared_ptr<MolVSValidations>> validations = {
      boost::make_shared<IsotopeValidation>()};
  MolVSValidation vm(validations);

  // testing MolVSDefault
  // testing for molecule with no atoms
  string smi1 = "";
  unique_ptr<ROMol> m1(SmilesToMol(smi1, 0, false));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (auto &query : errout1) {
    std::string msg = query.what();
    //    TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
    TEST_ASSERT(msg == "");
  }

  string smi2 = "O=C([O-])c1ccccc1";
  unique_ptr<ROMol> m2(SmilesToMol(smi2, 0, false));
  vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
  for (auto &query : errout2) {
    std::string msg = query.what();
    //    TEST_ASSERT(msg ==
    //                "INFO: [NeutralValidation] Not an overall neutral system
    //                (-1)");
    TEST_ASSERT(msg == "");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testAllowedAtomsValidation() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing AllowedAtoms validation"
      << std::endl;

  //	std::vector<string> atoms = {"C", "N", "O"};
  std::vector<unsigned int> atoms = {6, 7, 8};
  std::vector<shared_ptr<Atom>> atomList;

  for (auto &atom : atoms) {
    shared_ptr<Atom> a(new Atom(atom));
    atomList.push_back(a);
  }

  AllowedAtomsValidation vm(atomList);
  std::string smi1;

  smi1 = "CC(=O)CF";
  unique_ptr<ROMol> m1(SmilesToMol(smi1));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (auto &query : errout1) {
    std::string msg = query.what();
    TEST_ASSERT(
        msg ==
        "INFO: [AllowedAtomsValidation] Atom F is not in allowedAtoms list");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testDisallowedAtomsValidation() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing DisallowedAtoms validation"
      << std::endl;

  //	std::vector<string> atoms = {"F", "Cl", "Br"};
  std::vector<unsigned int> atoms = {9, 17, 35};
  std::vector<shared_ptr<Atom>> atomList;

  for (auto &atom : atoms) {
    shared_ptr<Atom> a(new Atom(atom));
    atomList.push_back(a);
  }

  DisallowedAtomsValidation vm(atomList);
  std::string smi1;

  smi1 = "CC(=O)CF";
  unique_ptr<ROMol> m1(SmilesToMol(smi1));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (auto &query : errout1) {
    std::string msg = query.what();
    TEST_ASSERT(
        msg ==
        "INFO: [DisallowedAtomsValidation] Atom F is in disallowedAtoms list");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testFragment() {
  BOOST_LOG(rdInfoLog)
      << "-----------------------\n Testing fragment validation" << std::endl;

  string smi1, smi2, smi3, smi4, smi5, smi6;
  MolVSValidation vm;

  // testing MolVSValidation fragmentValidation
  // FragmentValidation should identify 1,2-dichloroethane.
  smi1 = "ClCCCl.c1ccccc1O";
  unique_ptr<ROMol> m1(SmilesToMol(smi1, 0, false));
  vector<ValidationErrorInfo> errout1 = vm.validate(*m1, true);
  for (auto &query : errout1) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [FragmentValidation] 1,2-dichloroethane is present");
  }

  smi2 = "COCCOC.CCCBr";
  unique_ptr<ROMol> m2(SmilesToMol(smi2, 0, false));
  vector<ValidationErrorInfo> errout2 = vm.validate(*m2, true);
  for (auto &query : errout2) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [FragmentValidation] 1,2-dimethoxyethane is present");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testIs2DValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Is2DValidation"
                       << std::endl;

  Is2DValidation is2D;

  string mblock1 = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8753 4.9367 0 0
M  V30 2 C -0.4583 4.1667 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m1 {MolBlockToMol(mblock1, false, false)};
  vector<ValidationErrorInfo> errout1 = is2D.validate(*m1, true);
  TEST_ASSERT(errout1.empty());

  string mblock2 = R"(
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

  unique_ptr<ROMol> m2 {MolBlockToMol(mblock2, false, false)};
  vector<ValidationErrorInfo> errout2 = is2D.validate(*m2, true);
  TEST_ASSERT(errout2.size() == 1);
  string errmsg2 {errout2[0].what()};
  TEST_ASSERT(errmsg2 == "ERROR: [Is2DValidation] Molecule is 3D");

  string mblock3 = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 0 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m3 {MolBlockToMol(mblock3, false, false)};
  vector<ValidationErrorInfo> errout3 = is2D.validate(*m3, true);
  TEST_ASSERT(errout3.size() == 1);
  string errmsg3 {errout3[0].what()};
  TEST_ASSERT(errmsg3 == "ERROR: [Is2DValidation] All atoms have the same (x,y) coordinates");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testLayout2DValidation() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing Layout2DValidation"
                       << std::endl;
  Layout2DValidation layout2D;

  string mblock1 = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.8753 4.9367 0 0
M  V30 2 C -0.4583 4.1667 0 0
M  V30 3 C -0.2691 5.9671 0 0
M  V30 4 C -1.7337 5.4912 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 1 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m1 {MolBlockToMol(mblock1, false, false)};
  vector<ValidationErrorInfo> errout1 = layout2D.validate(*m1, true);
  TEST_ASSERT(errout1.empty());

  string mblock2 = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.6667 6.2067 0 0
M  V30 2 C -3.0004 5.4367 0 0
M  V30 3 C -3.0004 3.8965 0 0
M  V30 4 C -1.6667 3.1267 0 0
M  V30 5 C -0.3329 4.6000 0 0
M  V30 6 C -0.3329 4.7000 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 6
M  V30 3 1 2 3
M  V30 4 1 3 4
M  V30 5 1 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m2 {MolBlockToMol(mblock2, false, false)};
  vector<ValidationErrorInfo> errout2 = layout2D.validate(*m2, true);
  TEST_ASSERT(errout2.size() == 1);
  string errmsg2 {errout2[0].what()};
  TEST_ASSERT(errmsg2 == "ERROR: [Layout2DValidation] atom 5 too close to atom 6");

  string mblock3 = R"(
                    2D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.3333 4.8733 0 0
M  V30 2 C -2.5837 4.7283 0 0
M  V30 3 C -2.7087 3.4798 0 0
M  V30 4 C -1.6667 3.1267 0 0
M  V30 5 C -0.3329 3.8965 0 0
M  V30 6 C -0.9913 3.495 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 6
M  V30 3 1 2 3
M  V30 4 1 3 4
M  V30 5 1 4 5
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m3 {MolBlockToMol(mblock3, false, false)};
  vector<ValidationErrorInfo> errout3 = layout2D.validate(*m3, true);
  TEST_ASSERT(errout3.size() == 1);
  string errmsg3 {errout3[0].what()};
  TEST_ASSERT(errmsg3 == "ERROR: [Layout2DValidation] atom 6 too close to bond 5");

  string mblock4 = R"(
          01112413352D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -28.1663 10.4367 0 0
M  V30 2 C -29.5 9.6667 0 0
M  V30 3 C -29.5 11.2067 0 0
M  V30 4 F 0.0 10.4367 0 0
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

  unique_ptr<ROMol> m4 {MolBlockToMol(mblock4, false, false)};
  vector<ValidationErrorInfo> errout4 = layout2D.validate(*m4, true);
  TEST_ASSERT(errout4.size() == 1);
  string errmsg4 {errout4[0].what()};
  TEST_ASSERT(errmsg4 == "ERROR: [Layout2DValidation] length of bond 4 between atoms 1 and 4 exceeds a configured limit");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testValidateStereo() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing ValidateStereo"
                       << std::endl;

  StereoValidation stereo;
  string mblock, errmsg;
  vector<ValidationErrorInfo> errout;

  // 4 ligands - no issues
  mblock = R"(
          10052309532D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0.0003 7.27 0 0
M  V30 2 C -1.3333 6.5 0 0 CFG=1
M  V30 3 F -2.667 7.27 0 0
M  V30 4 O -1.3333 4.96 0 0
M  V30 5 C 0.0003 5.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 5 CFG=1
M  V30 2 1 2 3
M  V30 3 1 2 1
M  V30 4 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m0 {MolBlockToMol(mblock, false, false)};
  errout = stereo.validate(*m0, true);
  TEST_ASSERT(errout.size() == 0);

  // 4 ligands - too many stereo bonds with the same wedge/dash direction
  mblock = R"(
          10052310002D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0.0003 7.27 0 0
M  V30 2 C -1.3333 6.5 0 0 CFG=1
M  V30 3 F -2.667 7.27 0 0
M  V30 4 O -1.3333 4.96 0 0
M  V30 5 C 0.0003 5.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 5 CFG=1
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 1 CFG=1
M  V30 4 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m1 {MolBlockToMol(mblock, false, false)};
  Chirality::reapplyMolBlockWedging(*m1);
  errout = stereo.validate(*m1, true);
  TEST_ASSERT(errout.size() == 2);
  errmsg = errout[0].what();
  TEST_ASSERT(
    errmsg ==
    "ERROR: [StereoValidation] atom 2 has too many stereo bonds with like orientation");
  errmsg = errout[1].what();
  TEST_ASSERT(
    errmsg ==
    "ERROR: [StereoValidation] atom 2 has adjacent stereo bonds with like orientation");

  // 4 ligands - mismatching opposed wedge/dash bonds
  mblock = R"(
          10052311582D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0.0003 7.27 0 0
M  V30 2 C -1.3333 6.5 0 0
M  V30 3 F -2.667 7.27 0 0
M  V30 4 O -1.3333 4.96 0 0
M  V30 5 C 0.0003 5.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 5 CFG=1
M  V30 2 1 2 3 CFG=3
M  V30 3 1 2 1
M  V30 4 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  // 4 ligands - mismatching opposed wedge/dash bonds
  unique_ptr<ROMol> m2 {MolBlockToMol(mblock, false, false)};
  Chirality::reapplyMolBlockWedging(*m2);
  errout = stereo.validate(*m2, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0].what();
  TEST_ASSERT(
    errmsg ==
    "ERROR: [StereoValidation] atom 2 has opposing stereo bonds with different up/down orientation")

  // 4 ligands - potentially ambiguous umbrella configuration
  mblock = R"(
          10052313232D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0.0003 7.27 0 0
M  V30 2 C -1.3333 6.5 0 0
M  V30 3 F -2.667 7.27 0 0
M  V30 4 O -1.3333 4.96 0 0
M  V30 5 C 0.0003 5.73 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 5
M  V30 2 1 2 3 CFG=3
M  V30 3 1 2 1
M  V30 4 1 2 4
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m3 {MolBlockToMol(mblock, false, false)};
  Chirality::reapplyMolBlockWedging(*m3);
  errout = stereo.validate(*m3, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0].what();
  TEST_ASSERT(
    errmsg ==
    "ERROR: [StereoValidation] atom 2 has a potentially ambiguous representation: all non-stereo bonds opposite to the only stereo bond")

  // 4 ligands - colinearity / triangle rule violation
  mblock = R"(
          10052313312D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 F -1.083 6.5617 0 0
M  V30 2 C -2.4167 5.7917 0 0
M  V30 3 O -3.7503 6.5617 0 0
M  V30 4 C -1.2083 5.0017 0 0
M  V30 5 Cl -2.4167 6.5617 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 4 CFG=1
M  V30 2 1 2 1
M  V30 3 1 2 5
M  V30 4 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m4 {MolBlockToMol(mblock, false, false)};
  Chirality::reapplyMolBlockWedging(*m4);
  errout = stereo.validate(*m4, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0].what();
  TEST_ASSERT(
    errmsg ==
    "ERROR: [StereoValidation] colinearity or triangle rule violation of non-stereo bonds at atom 2")

  // 3 Ligands - No issues
  mblock = R"(
          10052313452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Cl -0.9997 6.895 0 0
M  V30 2 C -2.3333 6.125 0 0 CFG=2
M  V30 3 F -3.667 6.895 0 0
M  V30 4 C -2.3333 4.585 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 4 CFG=1
M  V30 2 1 2 3
M  V30 3 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m5 {MolBlockToMol(mblock, false, false)};
  Chirality::reapplyMolBlockWedging(*m5);
  errout = stereo.validate(*m5, true);
  TEST_ASSERT(errout.size() == 0);

  // 3 Ligands - multiple stereo bonds
  mblock = R"(
          10052313452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Cl -0.9997 6.895 0 0
M  V30 2 C -2.3333 6.125 0 0 CFG=2
M  V30 3 F -3.667 6.895 0 0
M  V30 4 C -2.3333 4.585 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 4 CFG=1
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m6 {MolBlockToMol(mblock, false, false)};
  Chirality::reapplyMolBlockWedging(*m6);
  errout = stereo.validate(*m6, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0].what();
  TEST_ASSERT(
    errmsg ==
    "ERROR: [StereoValidation] atom 2 has 3 explicit ligands and multiple stereo bonds")

  // 3 Ligands - colinearity violation
  mblock = R"(
          10052313452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Cl -0.9997 6.125 0 0
M  V30 2 C -2.3333 6.125 0 0 CFG=2
M  V30 3 F -3.667 6.125 0 0
M  V30 4 C -1.6667 4.71 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 4 CFG=1
M  V30 2 1 2 3
M  V30 3 1 2 1
M  V30 END BOND
M  V30 END CTAB
M  END
)";

  unique_ptr<ROMol> m7 {MolBlockToMol(mblock, false, false)};
  Chirality::reapplyMolBlockWedging(*m7);
  errout = stereo.validate(*m7, true);
  TEST_ASSERT(errout.size() == 1);
  errmsg = errout[0].what();
  TEST_ASSERT(
    errmsg ==
    "ERROR: [StereoValidation] colinearity of non-stereo bonds at atom 2");

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testValidateSmiles() {
  BOOST_LOG(rdInfoLog) << "-----------------------\n Testing ValidateSmiles"
                       << std::endl;

  // an invalid smiles should throw a ValueErrorException error
  try {
    vector<ValidationErrorInfo> errout1 = validateSmiles("3478q439g98h");
  } catch (const ValueErrorException &e) {
    std::string msg = e.what();
    TEST_ASSERT(msg ==
                "SMILES Parse Error: syntax error for input: 3478q439g98h")
  };

  vector<ValidationErrorInfo> errout2 = validateSmiles("");
  for (auto &query : errout2) {
    std::string msg = query.what();
    TEST_ASSERT(msg == "ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  vector<ValidationErrorInfo> errout3 = validateSmiles("ClCCCl.c1ccccc1O");
  for (auto &query : errout3) {
    std::string msg = query.what();
    TEST_ASSERT(msg ==
                "INFO: [FragmentValidation] 1,2-dichloroethane is present");
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main() {
  testRDKitValidation();
  testMolVSValidation();
  testMolVSOptions();
  testAllowedAtomsValidation();
  testDisallowedAtomsValidation();
  testFragment();
  testIs2DValidation();
  testLayout2DValidation();
  testValidateStereo();
  testValidateSmiles();
  return 0;
}
