// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ACTIVESUBSPACEREDUCER_HPP
#define ACTIVESUBSPACEREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"

namespace sgpp {
namespace base {

struct ActiveSubspaceInfo {
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
  std::vector<sgpp::base::DataMatrix> sampleMatrices;
};

class ActiveSubspaceReducer
    : public sgpp::base::FunctionReducer<ActiveSubspaceInfo> {
 public:
  class Gradient {
   public:
    virtual sgpp::base::DataVector gradientAt(sgpp::base::DataVector& v) = 0;
  };

  class GivenGradient : public Gradient {
   public:
    GivenGradient(std::shared_ptr<sgpp::optimization::VectorFunction> gradient);

    virtual sgpp::base::DataVector gradientAt(sgpp::base::DataVector& v);

   private:
    std::shared_ptr<sgpp::optimization::VectorFunction> gradient;
  };

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

  ActiveSubspaceReducer(size_t samples, std::shared_ptr<Gradient> gradient,
                        std::shared_ptr<VectorDistribution> distribution,
                        std::shared_ptr<CutoffCriterion<ActiveSubspaceInfo>> cutoff);

protected:
  void evaluateFunction(sgpp::optimization::ScalarFunction& input, ActiveSubspaceInfo& out);
  std::unique_ptr<sgpp::optimization::ScalarFunction> reduce(
      sgpp::optimization::ScalarFunction& input, size_t n, const ActiveSubspaceInfo& info) override;

 private:
  size_t samples;
  std::shared_ptr<Gradient> gradient;
  std::shared_ptr<VectorDistribution> distribution;
};

}  // namespace base
}  // namespace sgpp

#endif
