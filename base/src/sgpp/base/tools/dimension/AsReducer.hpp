// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASREDUCER_HPP
#define ASREDUCER_HPP

#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "EigenHelper.hpp"
#include "sgpp/base/function/vector/WrapperVectorFunction.hpp"

namespace sgpp {
namespace base {

struct AsInfo {
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
  DataMatrix permutation;
};

class TransfromFunction : public VectorFunction {
 public:
  TransfromFunction(const DataMatrix& m, size_t n) : VectorFunction(m.getNcols(), n), m(m) {
    this->m.resizeRowsCols(m.getNrows(), n);
    this->m.transpose();
  };
  ~TransfromFunction() override {};
  void eval(const DataVector& in, DataVector& out) override {
    out = in;
    out.sub(DataVector(in.size(), 0.5));
    if (out.max() != 0.0) {
      DataVector v = out;
      v.abs();
      out.mult(v.l2Norm() / v.sum());
    }
    out.mult(0.5);
    out = EigenHelper::mult(m, out);
    out.add(DataVector(out.size(), 0.5));
  };
  void clone(std::unique_ptr<VectorFunction>& clone) const override{};

 private:
  DataMatrix m;
};

template <class T>
struct AsResult : Result<T>{
  AsResult(const DataMatrix& m, size_t n) : f(m, n){
    DataMatrix transformation = m;
    transformation.transpose();
    transformation.resizeRowsCols(m.getNrows(), n);
    mInverse = transformation;
  }

  VectorFunction& getTransformationFunction() override
  { return f;
  };

 protected:

    void transformFrom(const DataVector& in, DataVector& out) {
    out = in;
    DataVector c(in.size(), 1);
      out.mult(2);
    out.sub(c);
    out = EigenHelper::mult(mInverse, out);
    out.mult(0.5);
    if (out.max() != 0.0) {
      DataVector v = out;
      v.abs();
      out.mult(v.sum() / 0.5);
    }
    out.add(DataVector(out.size(), 0.5));
  }

  DataMatrix mInverse;
  TransfromFunction f;
};

template <class I>
class AsReducer : public sgpp::base::Reducer<I, AsInfo, AsResult<I>> {
};

}  // namespace base
}  // namespace sgpp

#endif
