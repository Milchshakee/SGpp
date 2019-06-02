// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ACTIVESUBSPACEREDUCER_HPP
#define ACTIVESUBSPACEREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {

struct ActiveSubspaceInfo {
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
  Sample<sgpp::base::DataMatrix> sampleMatrices;
};

  struct ActiveSubspaceResult
  {
  DataMatrix transformation;

    std::unique_ptr<ReducedFunction> apply(
      sgpp::optimization::ScalarFunction& input);
  };

class ActiveSubspaceReducer
    : public sgpp::base::Reducer<Sample<DataMatrix>,ActiveSubspaceInfo, ActiveSubspaceResult> {
 public:

  static Sample<DataMatrix> fromGradientSample(const Sample<DataVector>& gradients);
  
  static Sample<DataMatrix> fromFiniteDifferences(optimization::ScalarFunction& func, VectorDistribution& v);

  class EigenValueCutoff : public CutoffCriterion<ActiveSubspaceInfo> {
   public:
    EigenValueCutoff(double minValue);

    size_t evaluate(const ActiveSubspaceInfo& info) override;

   private:
    double minEigenValue;
  };

  class IntervalCutoff : public CutoffCriterion<ActiveSubspaceInfo> {
   public:
    IntervalCutoff(size_t bootstrapSamples);

    size_t evaluate(const ActiveSubspaceInfo& info) override;

   private:
    size_t bootstrapSamples;
  };

  void evaluate(Sample<DataMatrix>& input, ActiveSubspaceInfo& out) override;
  ActiveSubspaceResult reduce(Sample<DataMatrix>& input, size_t c,
    const ActiveSubspaceInfo& info) override;
};

}  // namespace base
}  // namespace sgpp

#endif
