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
  Sample<AnovaBoundaryGrid::AnovaComponent, double> errorShares;
};

class AnovaResult : public Result<SGridSample> {
public:
  std::vector<bool> activeDimensions;
 size_t dimensions;
  double errorShare;

  AnovaResult() = default;
  AnovaResult(std::vector<bool>& activeDimensions, size_t dimensions, const SGridSample& sample,
              double errorShare);

  
    ScalarFunction& getOriginalFunction() override;
  VectorFunction& getTransformationFunction();
  SGridSample& getReducedOutput();
  ScalarFunction& getReducedFunctionSurrogate() override;

 private:
  EvalFunction originalFunction;
  EvalFunction eval;
  MatrixFunction f;
  SGridSample reducedSample;
};

class AnovaFixedCutter : public FixedCutter<SGridSample, AnovaInfo, AnovaResult> {
 public:
  AnovaFixedCutter(size_t n);

  AnovaResult cut(const SGridSample& input, const AnovaInfo& info) override;
};

  class AnovaErrorRuleCutter : public ErrorRuleCutter<SGridSample, AnovaInfo, AnovaResult> {
 public:
    AnovaErrorRuleCutter(ErrorRule& r, double maxError);

  AnovaResult cut(const SGridSample& input, const AnovaInfo& info) override;
  };

class AnovaReducer : public sgpp::base::Reducer<SGridSample, AnovaInfo, AnovaResult> {
 public:
  AnovaReducer(ErrorRule& rule);

  AnovaInfo evaluate(SGridSample& input) override;

private:
  ErrorRule& rule;
};

}  // namespace base
}  // namespace sgpp

#endif
