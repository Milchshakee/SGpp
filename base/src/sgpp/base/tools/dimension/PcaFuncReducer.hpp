// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/EvalFunction.hpp>
#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/base/tools/dimension/PcaReducer.hpp>

namespace sgpp {
namespace base {

struct PcaFuncInfo {
  sgpp::base::DataMatrix basis;
  DataVector mean;
};

struct PcaFuncInput {
  ScalarFunction& function;
  SGridSample& sample;
};

class PcaFuncResult : public Result<PcaFuncInput> {
 public:
  PcaFuncResult(const PcaFuncInput& input, const DataMatrix& m, size_t n, const DataVector& mean);

  ScalarFunction& getOriginalFunction() override;
  ScalarFunction& getReducedFunctionSurrogate() override;
  VectorFunction& getTransformationFunction() override;
  SGridSample& getReducedOutput() override;
  InputProjection& getProjection();

private:
  ScalarFunction* originalFunction;
  SGridSample reducedSurrogate;
 EvalFunction reducedSurrogateFunction;
  InputProjection projection;
};


class PcaFuncErrorRuleCutter : public ErrorRuleCutter<PcaFuncInput, PcaFuncInfo, PcaFuncResult> {
 public:
  PcaFuncErrorRuleCutter(ErrorRule& r, double maxError)
      : ErrorRuleCutter<PcaFuncInput, PcaFuncInfo, PcaFuncResult>(r, maxError) {
  }

  PcaFuncResult cut(const PcaFuncInput& input, const PcaFuncInfo& info) override;

};

class PcaFuncFixedCutter : public FixedCutter<PcaFuncInput, PcaFuncInfo, PcaFuncResult> {
 public:
  PcaFuncFixedCutter(size_t n);

  PcaFuncResult cut(const PcaFuncInput& input, const PcaFuncInfo& info) override;
};

class PcaFuncReducer : public Reducer<PcaFuncInput, PcaFuncInfo, PcaFuncResult> {
 public:
  PcaFuncReducer(std::shared_ptr<PcaSolver> solver, uint64_t seed, size_t samples, double stepSize,
                 size_t iterations);

  PcaFuncInfo evaluate(PcaFuncInput& input) override;

 private:
  std::shared_ptr<PcaSolver> solver;
  uint64_t seed;
  size_t samples;
  double stepSize;
  size_t iterations;

  void toDensity(SGridSample& input);
};

}  // namespace base
}  // namespace sgpp
