// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASMCREDUCER_HPP
#define ASMCREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "AsReducer.hpp"

namespace sgpp {
namespace base {

  struct AsMcInput
  {
  optimization::ScalarFunction& function;
  PointSample<DataMatrix> samples;
  };

  class AsMcResult : public AsResult<optimization::ScalarFunction>
  {
  public:
    AsMcResult(const AsMcInput& input, const DataMatrix& m, size_t n);


    AsMcResult(const AsMcResult& other)
      : AsResult<optimization::ScalarFunction>(other) {
      other.function->clone(function);
    }

    AsMcResult& operator=(const AsMcResult& other) {
      if (this == &other)
        return *this;
      AsResult<optimization::ScalarFunction>::operator=(other);
      other.function->clone(function);
      return *this;
    }

    optimization::ScalarFunction& getReducedFunction() override;
    optimization::ScalarFunction& getReducedOutput() override;

  private:
    std::unique_ptr<optimization::ScalarFunction> function;
  };

  class AsMcFixedCutter : public Cutter<AsMcInput, AsInfo, AsMcResult> {
 public:
  AsMcFixedCutter(size_t n);

  AsMcResult cut(const AsMcInput& input, const AsInfo& info) override;

 private:
  size_t n;
};

class AsMcIntervalCutter : public Cutter<AsMcInput, AsInfo, AsMcResult> {
 public:
  AsMcIntervalCutter(size_t bootstrapSamples);

  AsMcResult cut(const AsMcInput& input, const AsInfo& info) override;

 private:
  size_t bootstrapSamples;
};

class AsMcReducer : public sgpp::base::AsReducer<AsMcInput> {
 public:

  static PointSample<DataMatrix> fromGradientSample(const PointSample<DataVector>& gradients);
  
  static PointSample<DataMatrix> fromFiniteDifferences(optimization::ScalarFunction& func,
                                                       VectorDistribution& v);

  AsInfo evaluate(AsMcInput& input) override;
};

}  // namespace base
}  // namespace sgpp

#endif
