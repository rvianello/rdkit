import unittest

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import rdMinhash


class TestCase(unittest.TestCase):

  def testMinhashFromSparseFingerprint(self):
    mol1 = Chem.MolFromSmiles('O=C(O)CC1CCC1')
    mol2 = Chem.MolFromSmiles('O=C(O)CC1CC1')

    fpgen = rdFingerprintGenerator.GetMorganGenerator()
    bfp1, bfp2 = (fpgen.GetSparseFingerprint(mol) for mol in (mol1, mol2))

    # Hash1
    siggen = rdMinhash.MinhashSignatureGenerator32H1(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.7578125)

    siggen = rdMinhash.MinhashSignatureGenerator16H1(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.7578125)

    siggen = rdMinhash.MinhashSignatureGenerator8H1(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.765625)

    # Hash2
    siggen = rdMinhash.MinhashSignatureGenerator32H2(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.73046875)

    siggen = rdMinhash.MinhashSignatureGenerator16H2(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.73046875)

    siggen = rdMinhash.MinhashSignatureGenerator8H2(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.73828125)

  def testMinhashFromFingerprint(self):
    mol1 = Chem.MolFromSmiles('O=C(O)CC1CCC1')
    mol2 = Chem.MolFromSmiles('O=C(O)CC1CC1')

    fpgen = rdFingerprintGenerator.GetMorganGenerator(fpSize=2048)
    bfp1, bfp2 = (fpgen.GetFingerprint(mol) for mol in (mol1, mol2))

    # Hash1
    siggen = rdMinhash.MinhashSignatureGenerator32H1(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.70703125)

    siggen = rdMinhash.MinhashSignatureGenerator16H1(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.70703125)

    siggen = rdMinhash.MinhashSignatureGenerator8H1(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.71875)

    # Hash2
    siggen = rdMinhash.MinhashSignatureGenerator32H2(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.6875)

    siggen = rdMinhash.MinhashSignatureGenerator16H2(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.6875)

    siggen = rdMinhash.MinhashSignatureGenerator8H2(256)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))
    similarity = rdMinhash.TanimotoSimilarity(sig1, sig2)
    self.assertAlmostEqual(similarity, 0.7109375)

  def testLocalitySensitiveHashKeys(self):
    mol1 = Chem.MolFromSmiles('O=C(O)CC1CC1')
    mol2 = Chem.MolFromSmiles('O=C(O)CC1CCC1')

    fpgen = rdFingerprintGenerator.GetMorganGenerator(2, False)
    bfp1, bfp2 = (fpgen.GetSparseFingerprint(mol) for mol in (mol1, mol2))

    bands = 20
    rows = 5
    siggen = rdMinhash.MinhashSignatureGenerator32H2(bands*rows)
    sig1, sig2 = (siggen(bfp) for bfp in (bfp1, bfp2))

    keys1, keys2 = [
      rdMinhash.LocalitySensitiveHashKeys(bands, rows, sig)
      for sig in (sig1, sig2)
    ]
    self.assertEqual(len(keys1), bands)
    self.assertEqual(len(keys2), bands)

    matching = set(keys1) & set(keys2)
    self.assertEqual(len(matching), 7)


if __name__ == '__main__':

  unittest.main()
