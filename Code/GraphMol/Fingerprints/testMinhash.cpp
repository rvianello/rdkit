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
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/Minhash.h>
#include <DataStructs/BitOps.h>

#include <cmath>
#include <memory>

using namespace RDKit;

template <typename OutputType, typename HashType>
double signatureSimilarity(const SparseBitVect & fp1, const SparseBitVect & fp2, std::uint32_t l)
{
  using SignatureGenerator = Minhash::MinhashSignatureGenerator<OutputType, HashType>;
  using Signature = Minhash::MinhashSignature<OutputType>;

  SignatureGenerator signatureGenerator(l);
  Signature signature1 = signatureGenerator(fp1.getBitSet()->begin(), fp1.getBitSet()->end());
  Signature signature2 = signatureGenerator(fp2.getBitSet()->begin(), fp2.getBitSet()->end());

  return Minhash::TanimotoSimilarity(signature1, signature2);
}


void testBasic()
{
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Basic MinhashSignature Test" << std::endl;
  {
    auto pFpGenerator = std::unique_ptr<FingerprintGenerator<std::uint32_t>>(
      MorganFingerprint::getMorganGenerator<std::uint32_t>(2, false)
    );

    auto pMol1 = std::unique_ptr<ROMol>(SmilesToMol("O=C(O)CC1CC1"));
    auto pFp1 = std::unique_ptr<SparseBitVect>(pFpGenerator->getSparseFingerprint(*pMol1));
    TEST_ASSERT(pFp1->getNumOnBits() == 16);

    auto pMol2 = std::unique_ptr<ROMol>(SmilesToMol("O=C(O)CC1CCC1"));
    auto pFp2 = std::unique_ptr<SparseBitVect>(pFpGenerator->getSparseFingerprint(*pMol2));
    TEST_ASSERT(pFp2->getNumOnBits() == 18);

    BOOST_LOG(rdErrorLog)
      << "    fingerprint similarity: "
      << TanimotoSimilarity(*pFp1, *pFp2) << std::endl;

    auto sim32_h1 = signatureSimilarity<std::uint32_t, Minhash::Hash1>(*pFp1, *pFp2, 256);
    BOOST_LOG(rdErrorLog) << "    signature similarity (32, Hash1):" << sim32_h1 << std::endl;
    TEST_ASSERT(std::fabs(sim32_h1 - 0.8125) < 1e-4);

    auto sim16_h1 = signatureSimilarity<std::uint16_t, Minhash::Hash1>(*pFp1, *pFp2, 256);
    BOOST_LOG(rdErrorLog) << "    signature similarity (16, Hash1):" << sim16_h1 << std::endl;
    TEST_ASSERT(std::fabs(sim16_h1 - 0.8125) < 1e-4);
    
    auto sim8_h1 = signatureSimilarity<std::uint8_t, Minhash::Hash1>(*pFp1, *pFp2, 256);
    BOOST_LOG(rdErrorLog) << "    signature similarity ( 8, Hash1):" << sim8_h1 << std::endl;
    TEST_ASSERT(std::fabs(sim8_h1 - 0.816406) < 1e-6);
    
    auto sim32_h2 = signatureSimilarity<std::uint32_t, Minhash::Hash2>(*pFp1, *pFp2, 256);
    BOOST_LOG(rdErrorLog) << "    signature similarity (32, Hash2):" << sim32_h2 << std::endl;
    TEST_ASSERT(std::fabs(sim32_h2 - 0.800781) < 1e-6);
    
    auto sim16_h2 = signatureSimilarity<std::uint16_t, Minhash::Hash2>(*pFp1, *pFp2, 256);
    BOOST_LOG(rdErrorLog) << "    signature similarity (16, Hash2):" << sim16_h2 << std::endl;
    TEST_ASSERT(std::fabs(sim16_h2 - 0.800781) < 1e-6);
    
    auto sim8_h2 = signatureSimilarity<std::uint8_t, Minhash::Hash2>(*pFp1, *pFp2, 256);
    BOOST_LOG(rdErrorLog) << "    signature similarity ( 8, Hash2):" << sim8_h2 << std::endl;
    TEST_ASSERT(std::fabs(sim8_h2 - 0.804688) < 1e-6);
  }
}

int main(int, char *[])
{
  RDLog::InitLogs();

  testBasic();

  return 0;
}
