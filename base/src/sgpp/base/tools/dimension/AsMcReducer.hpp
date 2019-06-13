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

class AsMcIntervalCutter: public Cutter<PointSample<DataMatrix>,AsInfo, AsResult> {
 public:
  AsMcIntervalCutter(size_t bootstrapSamples);

  AsResult cut(const PointSample<DataMatrix>& input, const AsInfo& info) override;

 private:
  size_t bootstrapSamples;
};

class AsMcReducer : public sgpp::base::AsReducer<PointSample<DataMatrix>> {
 public:

  static PointSample<DataMatrix> fromGradientSample(const PointSample<DataVector>& gradients);
  
  static PointSample<DataMatrix> fromFiniteDifferences(optimization::ScalarFunction& func,
                                                       VectorDistribution& v);

  AsInfo evaluate(PointSample<DataMatrix>& input) override;
};

}  // namespace base
}  // namespace sgpp

#endif
