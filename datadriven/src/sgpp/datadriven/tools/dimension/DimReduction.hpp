// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DIMREDUCTION_HPP
#define DIMREDUCTION_HPP
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/tools/Sample.hpp>

namespace sgpp {
namespace base {

  struct ActiveSubspaceInfo {
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
};

  class InputProjectionFunction : public VectorFunction {
   public:
    InputProjectionFunction(const DataMatrix& basis, size_t reducedDims);
    ~InputProjectionFunction() override;
    void eval(const DataVector& x, DataVector& value) override;
    void clone(std::unique_ptr<VectorFunction>& clone) const override;

    static std::shared_ptr<InputProjectionFunction> identity(size_t dims);

   private:
    DataMatrix localBasis;
    DataMatrix activeSubspaceMat;
  };

  struct ReductionResult {
    ReductionResult(){};
    ReductionResult(
        std::shared_ptr<ScalarFunction>& func,
        std::shared_ptr<ScalarFunction> replacement);

    std::shared_ptr<ScalarFunction> originalFunction;
    std::shared_ptr<ScalarFunction> replacementFunction;
    std::shared_ptr<ScalarFunction> errorFunction;

    double mcL2Error(DistributionSample& sample);
  };

  struct GridReductionResult : public ReductionResult {
    GridReductionResult(){};
    GridReductionResult(std::shared_ptr<ScalarFunction>& func, std::shared_ptr<VectorFunction>& t,
                        std::shared_ptr<ScalarFunction>& sampleFunction, SGridSample& sample,
                        size_t gridPoints, double l2Error);

    SGridSample reducedSample;
    std::shared_ptr<VectorFunction> transformation;
    std::shared_ptr<ScalarFunction> reducedFunction;
    size_t gridPoints;
    double l2Error;
  };

  struct AsReductionResult : public ReductionResult {
    AsReductionResult(
        std::shared_ptr<ScalarFunction>& func,
        std::shared_ptr<ScalarFunction> replacement,
        std::vector<GridReductionResult>& results);

    std::vector<GridReductionResult> reductions;
  };

  namespace DimReduction
  {
  struct RegressionConfig {
    RegressionConfig(size_t reducedDimension) : reducedDimension(reducedDimension) {}

    size_t reducedDimension;
    size_t gridLevel = 3;
    std::vector<double> lambdas = {0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 1};
    std::vector<double> regularizationBases = {1.0, 0.5, 0.25, 0.125};
    size_t samples = 1000;
    size_t refinements = 1;
    size_t refinementPoints = 10;
    double trainDataShare = 0.02;
    size_t maxIterations = 1000;
    size_t crossValidations = 5;
    size_t subIterations = 3;
  };

  double calculateMcL2Error(ScalarFunction& func, DistributionSample& dist);

  ReductionResult reduceANOVA(ScalarFunction& f, const DataMatrix& basis, size_t reducedDims,
                              AnovaTypes::level_t level);
  sgpp::base::SGridSample createReducedAnovaSample(sgpp::base::SGridSample& sample,
                                                 AnovaTypes::level_t level, size_t reducedDims);

  AsReductionResult reduceAS(std::shared_ptr<ScalarFunction>& f,
                             sgpp::base::DistributionsVector dist,
                             double errorCalcSamples,
                           std::vector<RegressionConfig> config, bool useRelativeError = true);

  GridReductionResult activeSubspaceReductionStep(std::shared_ptr<ScalarFunction>& f,
                                              PointSample<double>& sample, const DataMatrix& basis,
                                              size_t reducedDims, RegressionConfig config,
                                              double totalError);

    PointSample<double> createActiveSubspaceSample(PointSample<double> input,
                                                   const DataMatrix& basis, size_t reducedDims);
  }
}  // namespace base
}  // namespace sgpp

#endif
