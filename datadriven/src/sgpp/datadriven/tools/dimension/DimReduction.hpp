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
    ReductionResult(
        std::shared_ptr<ScalarFunction>& func,
        std::shared_ptr<ScalarFunction> replacement);

    std::shared_ptr<ScalarFunction> originalFunction;
    std::shared_ptr<ScalarFunction> replacementFunction;
    std::shared_ptr<ScalarFunction> errorFunction;

    double mcL2Error(DistributionSample& sample);
  };

  struct GridReductionResult : public ReductionResult {
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
    double errorShareCutoff = 0.95;
    size_t gridLevel = 3;
    std::vector<double> lambdas = {0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 1};
    std::vector<double> regularizationBases = {1.0, 0.5, 0.25, 0.125};
    size_t samples = 1000;
    size_t refinements = 0;
    size_t refinementPoints = 0;
    double trainDataShare = 0.7;
    size_t maxIterations = 1000;
    size_t errorCalcSamples = 10000;
  };

double calculateMcL2Error(ScalarFunction& func,
                                                    VectorFunction& transformation,
                            ScalarFunction& reduced, DistributionSample& dist);
  double calculateMcL2Error(ScalarFunction& func, DistributionSample& dist);

  sgpp::base::SGridSample createReducedAnovaSample(sgpp::base::SGridSample& sample,
                                                 AnovaTypes::level_t level, size_t reducedDims);
  GridReductionResult activeSubspaceReduction(std::shared_ptr<ScalarFunction>& f, datadriven::Dataset& data, const DataMatrix& basis,
                                              size_t reducedDims, RegressionConfig config,
                                              double totalError);
  ReductionResult reduce(ScalarFunction& f, const DataMatrix& basis, size_t reducedDims,
                         AnovaTypes::level_t level);
  AsReductionResult reduceAS(std::shared_ptr<ScalarFunction>& f, sgpp::base::DistributionsVector dist,
                           std::vector<RegressionConfig> config);

    PointSample<double> createActiveSubspaceSample(PointSample<double> input,
                                                   const DataMatrix& basis, size_t reducedDims);
  }


template <class INPUT, class INFO, class OUTPUT>
class Cutter {
 public:
  virtual OUTPUT cut(const INPUT& input, const INFO& info) = 0;
};

class ErrorRule {
public:
  double calculateRelativeError(ScalarFunction& f, VectorFunction& t,
                                        ScalarFunction& r);
  virtual double calculateAbsoluteError(ScalarFunction& f, VectorFunction& t,
                                        ScalarFunction& r) = 0;
  virtual double calculateAbsoluteError(ScalarFunction& f) = 0;
};

    template <class INPUT, class INFO, class OUTPUT>
  class FixedCutter : public Cutter<INPUT, INFO, OUTPUT> {
   public:
    FixedCutter(size_t n)
      : n(n) {
    }

   protected:
    size_t n;
  };

template <class T>
class Result {
 public:
  double calculateRelativeError(ErrorRule& c)
  {
    return c.calculateRelativeError(getOriginalFunction(), getTransformationFunction(),
                            getReducedFunctionSurrogate());
  }

    double calculateAbsoluteError(ErrorRule& c) {
    return c.calculateAbsoluteError(getOriginalFunction(), getTransformationFunction(),
                            getReducedFunctionSurrogate());
  }

  virtual ScalarFunction& getOriginalFunction() = 0;
  virtual ScalarFunction& getReducedFunctionSurrogate() = 0;
  virtual VectorFunction& getTransformationFunction() = 0;
  virtual SGridSample& getReducedOutput() = 0;
};

  template <class INPUT, class INFO, class OUTPUT>
  class ErrorRuleCutter : public Cutter<INPUT, INFO, OUTPUT> {
 public:
  ErrorRuleCutter(ErrorRule& r, double maxError) : r(r), minVariance(maxError) {};

 protected:
  ErrorRule& r;
  double minVariance;
};

template <class INPUT, class INFO, class OUTPUT>
class Reducer {
 public:
  Reducer() = default;
  virtual ~Reducer() = default;

  OUTPUT evaluateAndCut(INPUT& input, Cutter<INPUT, INFO, OUTPUT>& cutter) {
    INFO i = evaluateFunction(input);
    OUTPUT result = cutter.cut(input, i);
    return result;
  }

  virtual INFO evaluate(INPUT& input) = 0;
};

  class InputProjection
  {
  public:
    InputProjection(const DataMatrix& basis, size_t n, const DataVector& mean);

    InputProjection& operator=(const InputProjection& other) {
      if (this == &other) return *this;
      oldToNewBasis = other.oldToNewBasis;
      newToOldBasis = other.newToOldBasis;
      oldDimensions = other.oldDimensions;
      newDimensions = other.newDimensions;
      mean = other.mean;
      posRange = other.posRange;
      negRange = other.negRange;
      start = other.start;
      end = other.end;
      func = ProjectionFunction(*this);
      return *this;
    }

    void inverse(const DataVector& in, DataVector& out);
    VectorFunction& getFunction();

    size_t getNewDimensions();
    const DataMatrix getTransformationMatrix();
    const DataVector& getStart();
    const DataVector& getEnd();

  private:
    DataMatrix oldToNewBasis;
    /// basis transposed and then cut
    DataMatrix newToOldBasis;
    size_t oldDimensions;
    size_t newDimensions;
    DataVector mean;
    DataVector posRange;
    DataVector negRange;
    DataVector start;
    DataVector end;

    void calculateRanges();

    class ProjectionFunction : public VectorFunction {
     public:
      ProjectionFunction() = default;
      ProjectionFunction(InputProjection& p);

      ~ProjectionFunction() = default;
      void eval(const DataVector& in, DataVector& out) override;
      void clone(std::unique_ptr<VectorFunction>& clone) const override;

     private:
      InputProjection* p;
    };

    ProjectionFunction func;
  };
}  // namespace base
}  // namespace sgpp

#endif
