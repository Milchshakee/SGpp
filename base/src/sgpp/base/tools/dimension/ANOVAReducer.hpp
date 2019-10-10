// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ANOVAREDUCER_HPP
#define ANOVAREDUCER_HPP

#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>
#include <sgpp/base/tools/Sample.hpp>
#include <sgpp/base/function/scalar/EvalFunction.hpp>
#include <sgpp/base/function/vector/MatrixFunction.hpp>

namespace sgpp {
namespace base {

struct AnovaInfo {
  Sample<AnovaBoundaryGrid::AnovaComponent, double> variances;
};

  struct AnovaInput
  {
  const SGridSample sample;
    ScalarFunction& originalFunction;
  };

class AnovaResult : public Result<SGridSample> {
public:
  std::vector<bool> activeDimensions;
 size_t dimensions;
  double coveredVariance;

  AnovaResult() = default;
  AnovaResult(std::vector<bool>& activeDimensions, size_t dimensions, const AnovaInput& input,
              double converedVariance);

  
    ScalarFunction& getOriginalFunction() override;
  VectorFunction& getTransformationFunction();
  SGridSample& getReducedOutput();
  ScalarFunction& getReducedFunctionSurrogate() override;

 private:
  ScalarFunction* originalFunction;
  EvalFunction eval;
  MatrixFunction transformation;
  SGridSample reducedSample;
};

class AnovaFixedCutter : public FixedCutter<AnovaInput, AnovaInfo, AnovaResult> {
 public:
  AnovaFixedCutter(size_t n);

  AnovaResult cut(const AnovaInput& input, const AnovaInfo& info) override;
};

  class AnovaErrorRuleCutter : public ErrorRuleCutter<AnovaInput, AnovaInfo, AnovaResult> {
 public:
    AnovaErrorRuleCutter(ErrorRule& r, double maxError);

  AnovaResult cut(const AnovaInput& input, const AnovaInfo& info) override;
  };

class AnovaReducer : public sgpp::base::Reducer<AnovaInput, AnovaInfo, AnovaResult> {
 public:
  AnovaReducer(ErrorRule& rule);

  AnovaInfo evaluate(AnovaInput& input) override;

private:
  ErrorRule& rule;
};

}  // namespace base
}  // namespace sgpp

#endif
