// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASREDUCER_HPP
#define ASREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "EigenHelper.hpp"

namespace sgpp {
namespace base {

struct AsInfo {
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
  DataMatrix permutation;
};

template <class T>
struct AsResult : Result<T>{
  AsResult(const DataMatrix& m, size_t n) {
    DataMatrix transformation = m;
    transformation.resizeRowsCols(m.getNrows(), n);
    mInverse = transformation;
    transformation.transpose();
    this->m = transformation;
    f = TransformationFunction(transformation);
  }

  TransformationFunction& getTransformationFunction() override
  { return f;
  };

 protected:
  void transformTo(DataVector& in, DataVector& out)
  { out = in;
    out.sub(DataVector(in.size(), 0.5));
    if (out.max() != 0.0) {
      out.mult(0.5 / out.max());
      }
    out = EigenHelper::mult(m, out);
      out.add(DataVector(out.size(), 0.5));
  }

    void transformFrom(DataVector& in, DataVector& out) {
    out = in;
    out = EigenHelper::mult(mInverse, out);
    if (out.max() != 0.0) {
      out.mult(1.0 / out.max());
    }
    out.add(DataVector(in.size(), 0.5));
  }

  DataMatrix m;
  DataMatrix mInverse;
  TransformationFunction f;
};

template <class I>
class AsReducer : public sgpp::base::Reducer<I, AsInfo, AsResult<I>> {
};

}  // namespace base
}  // namespace sgpp

#endif
