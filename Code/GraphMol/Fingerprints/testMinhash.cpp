//
//  2023, Riccardo Vianello
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/test.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/Minhash.h>
#include <DataStructs/BitOps.h>

//#include <GraphMol/FileParsers/MolSupplier.h>
//#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;

void testBasic()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test Basic MinhashSignature" << std::endl;
  {
    FingerprintGenerator<std::uint32_t> *radius2Generator =
        MorganFingerprint::getMorganGenerator<std::uint32_t>(2, false);

    ROMol *mol1 = SmilesToMol("O=C(O)CC1CC1");
    ExplicitBitVect *fp1 = radius2Generator->getFingerprint(*mol1);
    BOOST_LOG(rdErrorLog) << "    getNumBits: " << fp1->getNumBits() << std::endl;
    BOOST_LOG(rdErrorLog) << "    getNumOnBits: " << fp1->getNumOnBits() << std::endl;
    TEST_ASSERT(fp1->getNumOnBits() == 15);
    IntVect onBits1;
    fp1->getOnBits(onBits1);
    TEST_ASSERT(onBits1.size() == 15);

    ROMol *mol2 = SmilesToMol("O=C(OC)CC1CCC1");
    ExplicitBitVect *fp2 = radius2Generator->getFingerprint(*mol2);
    BOOST_LOG(rdErrorLog) << "    getNumOnBits: " << fp2->getNumOnBits() << std::endl;
    TEST_ASSERT(fp2->getNumOnBits() == 21);
    IntVect onBits2;
    fp2->getOnBits(onBits2);
    TEST_ASSERT(onBits2.size() == 21);

    BOOST_LOG(rdErrorLog)
      << "    fingerprint similarity: "
      << TanimotoSimilarity(*fp1, *fp2) << std::endl;

    using SignatureGenerator32 = Minhash::MinhashSignatureGenerator<uint32_t>;
    using Signature32 = SignatureGenerator32::MinhashSignature;

    SignatureGenerator32 minhashSignatureGenerator32(64);
    Signature32 signature32_1 = minhashSignatureGenerator32(onBits1.begin(), onBits1.end());
    Signature32 signature32_2 = minhashSignatureGenerator32(onBits2.begin(), onBits2.end());

    BOOST_LOG(rdErrorLog)
      << "    signature32 similarity: "
      << Minhash::similarity(signature32_1, signature32_2) << std::endl;
    //TEST_ASSERT(Minhash::similarity(signature32_1, signature32_2) == 0.414062);

    using SignatureGenerator8 = Minhash::MinhashSignatureGenerator<uint8_t>;
    using Signature8 = SignatureGenerator8::MinhashSignature;

    SignatureGenerator8 minhashSignatureGenerator8(256);
    Signature8 signature8_1 = minhashSignatureGenerator8(onBits1.begin(), onBits1.end());
    Signature8 signature8_2 = minhashSignatureGenerator8(onBits2.begin(), onBits2.end());

    BOOST_LOG(rdErrorLog)
      << "    signature8  similarity: "
      << Minhash::similarity(signature8_1, signature8_2) << std::endl;
    TEST_ASSERT(Minhash::similarity(signature8_1, signature8_2) == 0.359375);

    delete fp1;
    delete fp2;
    delete radius2Generator;
  }
}

int main(int, char *[])
{
  RDLog::InitLogs();

  testBasic();

  return 0;
}
