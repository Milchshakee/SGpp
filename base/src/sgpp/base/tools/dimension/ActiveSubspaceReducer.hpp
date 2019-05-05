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

class ActiveSubspaceReducer : public sgpp::base::FunctionReducer {
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

  class CutoffCriterion {
   public:
    virtual size_t evaluate(const sgpp::base::DataMatrix& eigenVectors,
                            const sgpp::base::DataVector& eigenValues,
                            const std::vector<sgpp::base::DataMatrix>& sampleMatrices) = 0;
  };

  class FixedCutoff : public CutoffCriterion {
   public:
    FixedCutoff(size_t n);

    size_t evaluate(const sgpp::base::DataMatrix& eigenVectors,
                    const sgpp::base::DataVector& eigenValues,
                    const std::vector<sgpp::base::DataMatrix>& sampleMatrices);

   private:
    size_t n;
  };

  class EigenValueCutoff : public CutoffCriterion {
   public:
    EigenValueCutoff(double minValue);

    size_t evaluate(const sgpp::base::DataMatrix& eigenVectors,
                    const sgpp::base::DataVector& eigenValues,
                    const std::vector<sgpp::base::DataMatrix>& sampleMatrices);

   private:
    double minEigenValue;
  };

    class IntervalCutoff : public CutoffCriterion {
   public:
      IntervalCutoff(size_t bootstrapSamples);

    size_t evaluate(const sgpp::base::DataMatrix& eigenVectors,
                    const sgpp::base::DataVector& eigenValues,
                    const std::vector<sgpp::base::DataMatrix>& sampleMatrices);

   private:
    size_t bootstrapSamples;
  };

  ActiveSubspaceReducer(size_t samples, std::shared_ptr < Gradient>  gradient,
                        std::shared_ptr < VectorDistribution>  distribution,
                        std::shared_ptr < CutoffCriterion>  cutoff);

  std::unique_ptr<sgpp::optimization::ScalarFunction> reduceFunction(
      sgpp::optimization::ScalarFunction& input) override;

 private:
  size_t samples;
  std::shared_ptr<Gradient> gradient;
  std::shared_ptr<VectorDistribution> distribution;
  std::shared_ptr<CutoffCriterion> cutoff;
};

}  // namespace base
}  // namespace sgpp

#endif
