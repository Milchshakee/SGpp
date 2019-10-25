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
  AsMcResult(const AsMcInput& input, const DataMatrix& m, size_t n, GridType type, level_t l,
             const DataVector& mean);

  AsMcResult& operator=(const AsMcResult& other) {
    if (this == &other) return *this;
    Result<SGridSample>::operator =(other);
    originalFunction = other.originalFunction;
    projection = other.projection;
    reduced = other.reduced;
    evalFunc = EvalFunction(reduced);
    return *this;
  }

  ScalarFunction& getOriginalFunction() override;
  VectorFunction& getTransformationFunction() override;
  ScalarFunction& getReducedFunctionSurrogate() override;
  SGridSample& getReducedOutput() override;
  InputProjection& getProjection();

 private:
  ScalarFunction* originalFunction;
  InputProjection projection;
  SGridSample reduced;
  EvalFunction evalFunc;
};

class AsMcCutter {
public:
  AsMcCutter(GridType type, level_t level, const DataVector& mean);

 protected:
  GridType type;
  level_t level;
  DataVector mean;
};

class AsMcFixedCutter : public FixedCutter<AsMcInput, AsInfo, AsMcResult>, public AsMcCutter {
 public:
  AsMcFixedCutter(size_t n, GridType type, level_t level, const DataVector& mean);

  AsMcResult cut(const AsMcInput& input, const AsInfo& info) override;
};

class AsMcErrorRuleCutter : public ErrorRuleCutter<AsMcInput, AsInfo, AsMcResult>,
                            public AsMcCutter {
 public:
  AsMcErrorRuleCutter(ErrorRule& r, double maxError, GridType type, level_t level,
                      const DataVector& mean);

  AsMcResult cut(const AsMcInput& input, const AsInfo& info) override;
};

class AsMcIntervalCutter : public Cutter<AsMcInput, AsInfo, AsMcResult>, public AsMcCutter {
 public:
  AsMcIntervalCutter(size_t bootstrapSamples, GridType type, level_t level, const DataVector& mean);

  AsMcResult cut(const AsMcInput& input, const AsInfo& info) override;

 private:
  size_t bootstrapSamples;
};

class AsMcReducer : public Reducer<AsMcInput, AsInfo, AsMcResult> {
 public:
  static PointSample<DataMatrix> fromGradientSample(const PointSample<DataVector>& gradients);

  static PointSample<DataMatrix> fromFiniteDifferences(ScalarFunction& func, VectorDistribution& v,
                                                       double h);

  AsInfo evaluate(AsMcInput& input) override;
};

}  // namespace base
}  // namespace sgpp

#endif
