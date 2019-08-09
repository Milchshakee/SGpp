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
};

class AnovaResult : public Result<SGridSample> {
public:
  std::vector<bool> activeDimensions;
 size_t dimensions;

  AnovaResult() = default;
  AnovaResult(std::vector<bool>& activeDimensions,
             size_t dimensions, const SGridSample& sample);


  ScalarFunction& getOriginalFunction() override;
  VectorFunction& getTransformationFunction();
  SGridSample& getReducedOutput();
  ScalarFunction& getReducedFunction() override;

 private:
  EvalFunction originalFunction;
  EvalFunction eval;
  MatrixFunction f;
  SGridSample reducedSample;
};

class AnovaFixedCutter : public FixedCutter<SGridSample, AnovaInfo, AnovaResult> {
 public:
  AnovaFixedCutter(ErrorRule& r, size_t n);

  AnovaResult cut(const SGridSample& input, const AnovaInfo& info) override;
};

  class AnovaErrorRuleCutter : public ErrorRuleCutter<SGridSample, AnovaInfo, AnovaResult> {
 public:
    AnovaErrorRuleCutter(ErrorRule& r, double maxError);

  AnovaResult cut(const SGridSample& input, const AnovaInfo& info) override;
  };

class AnovaReducer : public sgpp::base::Reducer<SGridSample, AnovaInfo, AnovaResult> {
 public:
  AnovaInfo evaluate(SGridSample& input) override;
};

}  // namespace base
}  // namespace sgpp

#endif
