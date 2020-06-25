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
        std::shared_ptr<ScalarFunction> replacement);

    std::shared_ptr<ScalarFunction> replacementFunction;
  };

  struct GridReductionResult : public ReductionResult {
    GridReductionResult(){};
    GridReductionResult(std::shared_ptr<VectorFunction>& t,
                        std::shared_ptr<ScalarFunction>& sampleFunction, SGridSample& sample);

    SGridSample reducedSample;
    std::shared_ptr<VectorFunction> transformation;
    std::shared_ptr<ScalarFunction> reducedFunction;
  };

  struct AsReductionResult : public ReductionResult {
    AsReductionResult(
        std::shared_ptr<ScalarFunction> replacement,
        std::vector<GridReductionResult>& results);

    std::vector<GridReductionResult> reductions;
  };

  namespace DimReduction
  {
  struct ReductionConfig {
    ReductionConfig(size_t maxReducedDimension) : maxReducedDimension(maxReducedDimension) {}

    size_t maxReducedDimension;
    size_t maxIterations = 10;
    double maxErrorShare = 0.01;
    size_t errorCalcSamples = 10000;
  };

  struct RegressionConfig {

    size_t maxGridPoints = 100;
    std::vector<double> lambdas = {0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 1};
    std::vector<double> regularizationBases = {1.0, 0.5, 0.25, 0.125};
    size_t refinements = 1;
    double refinementPointsPerGridPoint = 0.1;
    double trainDataPerGridPoint = 1;
    size_t crossValidations = 5;
    size_t sampleCount = 10000;
  };

    struct GradientConfig
    {
    enum Type {
      GRADIENT_FUNCTION,
      NEAREST_NEIGHBOUR_APPROXIMATION,
      RANDOM_NEIGHBOUR_APPROXIMATION,
    };
    Type type;

    size_t sampleCount = 10000;

    std::shared_ptr<VectorFunction> initialGradientFunction;
    double h = 1e-8;

    double pNorm = 1;
    double maxNorm = 0.2;

    size_t neighbourConsiderationCount = 1000;
    };

  struct BasisConfig {
    enum Type {
      ACTIVE_SUBSPACE,
      INV_AS,
      RANDOM
    };
    std::vector<Type> types;

    double eigenValueShare = 0.9;
    bool useBiggestInterval = true;
    size_t basisIterations = 10;
  };

  struct Error {
    Error(std::string name) : name(name) {}

    virtual double calculateError(std::shared_ptr<ScalarFunction>& func,
                                  PointSample<double>& errorSample) = 0;

    const std::string& getName()
    { return name;
    }

  private:
    std::string name;
  };

    struct ExaminationConfig
    {
    ExaminationConfig(Error& e) : error(e) {}

    bool hasOutput()
    { return !outputDir.empty();
    }

      Error& error;
      bool useRelativeError = true;
      uint64_t randomSeed;
    bool discardWorseIterations = false;
      std::string outputDir;
    };

    struct Output
    {
      DataMatrix eigenVectors;
      DataVector eigenValues;
      size_t gridPointsAdded;
      double newError;
    };

    class L2 : public Error
    {
    public:
      L2()
        : Error("RMSE") {
      }

      double calculateError(std::shared_ptr<ScalarFunction>& func,
        PointSample<double>& errorSample) override;
    };

    class MSE : public Error
    {
     public:
      MSE() : Error("MSE") {}
      double calculateError(std::shared_ptr<ScalarFunction>& func,
                            PointSample<double>& errorSample) override;
    };

    class NRMSE : public Error {
     public:
      NRMSE() : Error("NRMSE") {}
      double calculateError(std::shared_ptr<ScalarFunction>& func,
                            PointSample<double>& errorSample) override;
    };


    extern  L2 L2_ERROR;
    extern L2 RMSE_ERROR;
    extern MSE MSE_ERROR;
    extern  NRMSE NRMSE_ERROR;

  std::shared_ptr<VectorFunction> finiteDifferencesFunction(std::shared_ptr<ScalarFunction>& f,
                                                            double h);

  AsReductionResult reduceAS(PointSample<double>& sample, ReductionConfig& redConfig,
                             RegressionConfig& fastRegConfig,
                               RegressionConfig& regConfig, GradientConfig& gradConfig,
                               BasisConfig& basisConfig, ExaminationConfig& examConfig, std::vector<Output>& output);

  }
}  // namespace base
}  // namespace sgpp

#endif
