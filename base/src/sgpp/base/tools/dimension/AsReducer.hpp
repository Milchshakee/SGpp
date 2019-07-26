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

namespace sgpp {
namespace base {

struct AsInfo {
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
};

template <class T>
struct AsResult : Result<T>{
  AsResult(const DataMatrix& m, size_t n) {
    DataMatrix transformation = m;
    //transformation.resizeRowsCols(m.getNrows(), n);
    transformation.setColumn(1, DataVector(2, 0));
    transformation.transpose();
    f = TransformationFunction(transformation);
  }

  TransformationFunction& getTransformationFunction() override
  { return f;
  };

 private:
  TransformationFunction f;
};

template <class I>
class AsReducer : public sgpp::base::Reducer<I, AsInfo, AsResult<I>> {
};

}  // namespace base
}  // namespace sgpp

#endif
