//
//  2023, Riccardo Vianello
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/Wrap.h>
#include <boost/python.hpp>
#include <GraphMol/Fingerprints/Minhash.h>
#include <DataStructs/SparseBitVect.h>
#include <DataStructs/ExplicitBitVect.h>

#include <string>

namespace python = boost::python;

namespace RDKit {
namespace Minhash {

namespace {
  template <typename OutputType>
  double similarity(const MinhashSignature<OutputType> *s1, const MinhashSignature<OutputType> *s2) {
    return TanimotoSimilarity(*s1, *s2);
  }

  template <typename OutputType, typename HashType>
  MinhashSignature<OutputType> *
  callSparseBitVect(const MinhashSignatureGenerator<OutputType, HashType> *gen, const SparseBitVect *fp) {
    MinhashSignature<OutputType> signature = (*gen)(fp->getBitSet()->begin(), fp->getBitSet()->end());
    return new MinhashSignature<OutputType>(std::move(signature));
  }

  template <typename OutputType, typename HashType>
  MinhashSignature<OutputType> *
  callExplicitBitVect(const MinhashSignatureGenerator<OutputType, HashType> *gen, const ExplicitBitVect *fp) {
    IntVect v;
    fp->getOnBits(v);
    MinhashSignature<OutputType> signature = (*gen)(v.begin(), v.end());
    return new MinhashSignature<OutputType>(std::move(signature));
  }

  template <typename OutputType>
  OutputType signatureGetItem(const MinhashSignature<OutputType> & sig, int pos) {
    int size = static_cast<int>(sig.size());
    if (pos < 0) {
      pos += size;
      if (pos < 0) {
        throw IndexErrorException(pos);
      }
    }
    if (pos >= size) {
      throw IndexErrorException(pos);
    }
    return sig[pos];
  }

  template <typename OutputType>
  void wrapSignature(const std::string & nm) {
    python::class_<MinhashSignature<OutputType>>(nm.c_str(), python::no_init)
      .def("__len__", &MinhashSignature<OutputType>::size)
      .def("__getitem__", signatureGetItem<OutputType>, python::arg("pos"))
      ;
  }

  template <typename OutputType, typename HashType>
  void wrapGenerator(const std::string & nm) {
    python::class_<MinhashSignatureGenerator<OutputType, HashType>>(
      nm.c_str(), python::init<std::uint32_t, python::optional<int>>()
      )
      .def("__call__", callSparseBitVect<OutputType, HashType>,
           (python::arg("bfp"), python::return_value_policy<python::manage_new_object>()))
      .def("__call__", callExplicitBitVect<OutputType, HashType>,
           (python::arg("bfp"), python::return_value_policy<python::manage_new_object>()))
      ;
  }
  
}

BOOST_PYTHON_MODULE(rdMinhash) {
  wrapSignature<std::uint32_t>("MinhashSignature32");
  wrapSignature<std::uint16_t>("MinhashSignature16");
  wrapSignature<std::uint8_t>("MinhashSignature8");

  wrapGenerator<std::uint32_t, Hash1>("MinhashSignatureGenerator32H1");
  wrapGenerator<std::uint16_t, Hash1>("MinhashSignatureGenerator16H1");
  wrapGenerator<std::uint8_t, Hash1>("MinhashSignatureGenerator8H1");

  wrapGenerator<std::uint32_t, Hash2>("MinhashSignatureGenerator32H2");
  wrapGenerator<std::uint16_t, Hash2>("MinhashSignatureGenerator16H2");
  wrapGenerator<std::uint8_t, Hash2>("MinhashSignatureGenerator8H2");

  python::def("TanimotoSimilarity", similarity<std::uint32_t>, (python::arg("s1"), python::arg("s2")));
  python::def("TanimotoSimilarity", similarity<std::uint16_t>, (python::arg("s1"), python::arg("s2")));
  python::def("TanimotoSimilarity", similarity<std::uint8_t>, (python::arg("s1"), python::arg("s2")));
}

} // namespace Minhash
} // namespace RDKit
