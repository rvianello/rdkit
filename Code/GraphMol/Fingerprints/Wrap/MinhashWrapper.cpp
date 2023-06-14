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

#include <string>

namespace python = boost::python;

namespace RDKit {
namespace Minhash {

namespace {
  template <typename OutputType, typename HashType>
  typename MinhashSignatureGenerator<OutputType, HashType>::MinhashSignature
  call(const MinhashSignatureGenerator<OutputType, HashType> *gen, const SparseBitVect *fp) {
    return (*gen)(fp->getBitSet()->begin(), fp->getBitSet()->end());
  }

  template <typename OutputType, typename HashType>
  void wrapGenerator(const std::string & nm) {
    python::class_<MinhashSignatureGenerator<OutputType, HashType>>(
      nm.c_str(), python::init<std::uint32_t, python::optional<int>>()
      )
      .def("__call__", call<OutputType, HashType>, python::arg("bfp"))
      ;
  }
  
}

BOOST_PYTHON_MODULE(rdMinhash) {
  wrapGenerator<std::uint32_t, Hash1>("MinhashSignatureGenerator32H1");
  wrapGenerator<std::uint16_t, Hash1>("MinhashSignatureGenerator16H1");
  wrapGenerator<std::uint8_t, Hash1>("MinhashSignatureGenerator8H1");

  wrapGenerator<std::uint32_t, Hash2>("MinhashSignatureGenerator32H2");
  wrapGenerator<std::uint16_t, Hash2>("MinhashSignatureGenerator16H2");
  wrapGenerator<std::uint8_t, Hash2>("MinhashSignatureGenerator8H2");
}

} // namespace Minhash
} // namespace RDKit
