//
//  2023, Riccardo Vianello
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file Minhash.h

*/
#include <RDGeneral/export.h>
#ifndef RD_MINHASH_H_2023_06
#define RD_MINHASH_H_2023_06

#include <algorithm>
#include <cstdint>
#include <functional>
#include <map>
#include <random>
#include <set>
#include <vector>

namespace RDKit {
namespace Minhash {

// Families of universal hash functions
// https://en.wikipedia.org/wiki/Universal_hashing

inline std::uint64_t hash(std::uint64_t x, std::uint8_t l, std::uint64_t a) {
  // hashes x 2/m-almost-universally into l <= 64 bits
  // using random odd seed a
  return (a*x) >> (64-l);
};

inline std::uint32_t hash(std::uint32_t x, std::uint8_t l, std::uint64_t a, std::uint64_t b) {
  // hashes x truly universally into l <= 32 bits
  // using random seeds a and b
  // 0 < a < 2**64
  // 0 <= b < 2**64
  return (a*x+b) >> (64-l);
};

using RandomFunction = std::function<std::uint64_t ()>;

struct Hash1 {
  using InputType = std::uint64_t;
  using OutputType = std::uint64_t;
  using HashSeed = std::uint64_t;

  static HashSeed generate_seed(RandomFunction random) 
  {
    return random() | 1;
  }

  explicit Hash1(std::uint8_t l) : l{l} {}
  OutputType operator()(InputType x, HashSeed seed) {
    return hash(x, l, seed);
  }
  std::uint8_t l;
};

struct Hash2 {
  using InputType = std::uint32_t;
  using OutputType = std::uint32_t;
  using HashSeed = std::pair<std::uint64_t, std::uint64_t>;

  static HashSeed generate_seed(RandomFunction random) 
  {
    std::uint64_t b{random()};
    std::uint64_t a{random()};
    while (a == 0) {a = random();}
    return {a, b};
  }

  explicit Hash2(std::uint8_t l) : l{l} {}
  OutputType operator()(InputType x, HashSeed seed) {
    return hash(x, l, seed.first, seed.second);
  }
  std::uint8_t l;
};

template <typename T> class MinhashSignature : public std::vector<T> {};

template <typename OutputType, typename HashType>
class MinhashSignatureGenerator {
public:
  using SignatureType = MinhashSignature<OutputType>;

  MinhashSignatureGenerator(std::uint32_t size, int rng_seed=42)
    : seeds{} {
    PRECONDITION(size, "The signature size must be larger than zero");
  
    std::mt19937_64 rng(rng_seed);
    std::uniform_int_distribution<std::uint64_t> distribution;
    auto random = [&rng, &distribution]() {return distribution(rng);};

    // generate a collection of "size" unique seed parameters for the
    // hash functions
    std::set<typename HashType::HashSeed> temp;
    while (temp.size() < size) {
      temp.insert(HashType::generate_seed(random));
    }
    seeds.reserve(size);
    std::copy(temp.begin(), temp.end(), std::back_inserter(seeds));
  }

  template <typename InputIt>
  SignatureType operator()(InputIt beginIt, InputIt endIt) const
  {
    HashType hash{l};
    SignatureType signature;
    signature.reserve(seeds.size());
    // loop over the hash functions, for each function compute a minhash
    // append the minhash to the signature
    for (const auto & seed : seeds) {
      InputIt it{beginIt};
      OutputType minhash = hash(*it, seed);
      while (++it != endIt) {
        minhash = std::min(minhash, static_cast<OutputType>(hash(*it, seed)));
      }
      signature.push_back(minhash);
    }
    return signature;
  }

private:
  using HashSeeds = std::vector<typename HashType::HashSeed>;
  HashSeeds seeds;
  static constexpr std::uint8_t l = sizeof(OutputType)*8;
};

template <typename T>
double TanimotoSimilarity(const T & sign1, const T & sign2)
{
  PRECONDITION(sign1.size() == sign2.size(), "The input signatures must have the same length");
  PRECONDITION(sign1.size(), "The input signatures can't be empty");

  int same{};

  for (auto it1 = sign1.begin(), it2 = sign2.begin();
       it1 != sign1.end() && it2 != sign2.end();
       ++it1, ++it2
      )
  {
    if (*it1 == *it2) {
      same += 1;
    }
  }

  return static_cast<double>(same)/sign1.size();
}

template <typename T>
std::vector<std::uint32_t> LocalitySensitiveHashKeys(int bands, int rows, const T & signature)
{
  PRECONDITION(bands > 0, "The number of bands must > 0");
  PRECONDITION(rows > 0, "The number of rows must > 0");
  PRECONDITION(
    int(signature.size()) == bands*rows, "The length of the input signature must equal bands*rows");

  std::vector<std::uint32_t> keys;
  keys.reserve(bands);

  /*
   * compute the hash keys.
   * the rows from each band are first hashed together and then
   * prefixed with the band number to make sure that only the keys
   * from corresponding bands can match.
   */
  auto it = signature.begin();
  for (int band=0; band < bands; ++band) {
    // https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector/72073933#72073933
    uint32_t key = band * 53;
    for (int row=0; row < rows; ++row) {
      uint32_t x = *it++;
      x = ((x >> 16) ^ x) * 0x45d9f3b;
      x = ((x >> 16) ^ x) * 0x45d9f3b;
      x = (x >> 16) ^ x;
      key ^= x + 0x9e3779b9 + (key << 6) + (key >> 2);
    }
    // fold to 24 bits (http://isthe.com/chongo/tech/comp/fnv/#xor-fold)
    key = (key >> 24) ^ (key & (uint32_t)0xffffff);
    // put the band number in the most significant byte
    key |= (uint32_t)band << 24;

    keys.push_back(key);
  }

  return keys;
}

} // namespace Minhash
} // namespace RDKit

#endif