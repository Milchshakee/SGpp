// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ACTIVESUBSPACEREDUCER_HPP
#define ACTIVESUBSPACEREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "AsReducer.hpp"

namespace sgpp {
namespace base {

class AsMcIntervalCutter: public Cutter<AsInfo> {
 public:
  AsMcIntervalCutter(size_t bootstrapSamples);

  size_t evaluate(const AsInfo& info) override;

 private:
  size_t bootstrapSamples;
};

class AsMcReducer
    : public sgpp::base::AsReducer<Sample<DataMatrix>> {
 public:

  static Sample<DataMatrix> fromGradientSample(const Sample<DataVector>& gradients);
  
  static Sample<DataMatrix> fromFiniteDifferences(optimization::ScalarFunction& func, VectorDistribution& v);

  void evaluate(Sample<DataMatrix>& input, AsInfo& out) override;
  AsResult reduce(Sample<DataMatrix>& input, size_t c,
    const AsInfo& info) override;
};

}  // namespace base
}  // namespace sgpp

#endif
