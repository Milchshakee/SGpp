// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/tools/VectorDistribution.hpp>
#include "DimReduction.hpp"
#include "PcaReducer.hpp"
#include <sgpp/base/function/scalar/EvalFunction.hpp>

namespace sgpp {
namespace base {

struct PcaFuncInfo {
  sgpp::base::DataMatrix basis;
  sgpp::base::DataVector eigenValues;
  DataVector mean;
  DataVector varianceShares;
};

class PcaFuncResult : public Result<SGridSample> {
public:
  PcaFuncResult(const SGridSample& input, const DataMatrix& m, size_t n, const DataVector& mean);

  ScalarFunction& getReducedFunction() override;
  VectorFunction& getTransformationFunction() override;
  SGridSample& getReducedOutput() override;

private:
  SGridSample reduced;
 EvalFunction evalFunc;
  DataMatrix m;
  DataMatrix mInv;
 DataVector mean;

  void transformFrom(const DataVector& in, DataVector& out);
};

class PcaFuncCutter : public Cutter<SGridSample, PcaFuncInfo, PcaFuncResult> {};

class PcaFuncVarianceCutter : public PcaFuncCutter {
 public:
  PcaFuncVarianceCutter(double minVarianceShare);

  PcaFuncResult cut(const SGridSample& input, const PcaFuncInfo& info) override;

 private:
  double minVarianceShare;
};

class PcaFuncFixedCutter : public PcaFuncCutter {
 public:
  PcaFuncFixedCutter(size_t n);

  PcaFuncResult cut(const SGridSample& input, const PcaFuncInfo& info) override;

 private:
  size_t n;
};

class PcaFuncReducer : public Reducer<SGridSample, PcaFuncInfo, PcaResult> {
 public:
  PcaFuncReducer(std::shared_ptr<PcaSolver> solver, uint64_t seed, size_t samples, double stepSize,
                 size_t iterations);

  PcaFuncInfo evaluate(SGridSample& input) override;

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
