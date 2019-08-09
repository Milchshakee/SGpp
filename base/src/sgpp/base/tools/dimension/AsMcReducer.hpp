// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASMCREDUCER_HPP
#define ASMCREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/EvalFunction.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/tools/Sample.hpp>
#include <sgpp/base/tools/dimension/AsReducer.hpp>

namespace sgpp {
namespace base {

struct AsMcInput {
  ScalarFunction& function;
  PointSample<DataMatrix> samples;
};

class AsMcResult : public Result<SGridSample> {
 public:
  AsMcResult(const AsMcInput& input, const DataMatrix& m, size_t n, GridType type, level_t l);

  ScalarFunction& getOriginalFunction() override;
  VectorFunction& getTransformationFunction() override;
  ScalarFunction& getReducedFunction() override;
  SGridSample& getReducedOutput() override;

 private:
  InputProjection projection;
  SGridSample reduced;
  EvalFunction evalFunc;
  ScalarFunction* originalFunc;
};

class AsMcCutter {
public:
  AsMcCutter(GridType type, level_t level);

 protected:
  GridType type;
  level_t level;
};

class AsMcFixedCutter : public FixedCutter<AsMcInput, AsInfo, AsMcResult>, public AsMcCutter {
 public:
  AsMcFixedCutter(ErrorRule& r, size_t n, GridType type, level_t level);

  AsMcResult cut(const AsMcInput& input, const AsInfo& info) override;
};

class AsMcErrorRuleCutter : public ErrorRuleCutter<AsMcInput, AsInfo, AsMcResult>,
                            public AsMcCutter {
 public:
  AsMcErrorRuleCutter(ErrorRule& r, double maxError, GridType type, level_t level);

  AsMcResult cut(const AsMcInput& input, const AsInfo& info) override;
};

class AsMcIntervalCutter : public Cutter<AsMcInput, AsInfo, AsMcResult>, public AsMcCutter {
 public:
  AsMcIntervalCutter(size_t bootstrapSamples, GridType type, level_t level);

  AsMcResult cut(const AsMcInput& input, const AsInfo& info) override;

 private:
  size_t bootstrapSamples;
};

class AsMcReducer : public Reducer<AsMcInput, AsInfo, AsMcResult> {
 public:
  static PointSample<DataMatrix> fromGradientSample(const PointSample<DataVector>& gradients);

  static PointSample<DataMatrix> fromFiniteDifferences(ScalarFunction& func, VectorDistribution& v);

  AsInfo evaluate(AsMcInput& input) override;
};

}  // namespace base
}  // namespace sgpp

#endif
