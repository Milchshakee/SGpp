// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ANOVAREDUCER_HPP
#define ANOVAREDUCER_HPP

#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include "sgpp/base/grid/type/AnovaBoundaryGrid.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {

struct AnovaInfo {
  sgpp::base::Sample<AnovaHelper::AnovaComponent, double> variances;
  double totalVariance;
};

class AnovaResult : public Result<SGridSample> {
public:
  std::vector<bool> activeDimensions;
 AnovaHelper::AnovaComponentVector activeComponents;
  double coveredVariance;
 size_t dimensions;

  AnovaResult(std::vector<bool>& activeDimensions,
             AnovaHelper::AnovaComponentVector& activeComponents, double coveredVariance,
             size_t dimensions, SGridSample sample);

 double calcMcL2Error(optimization::ScalarFunction& func, size_t paths, uint64_t seed = std::mt19937_64::default_seed);

  TransformationFunction& getTransformationFunction();
  SGridSample& getReducedOutput();

 private:
  TransformationFunction f;
  SGridSample orginialSample;
  SGridSample reducedSample;
};

class AnovaCutter : public sgpp::base::Cutter<SGridSample, AnovaInfo, AnovaResult> {};

class AnovaFixedCutter : public AnovaCutter {
 public:
  AnovaFixedCutter(size_t n);

  AnovaResult cut(const SGridSample& input, const AnovaInfo& info) override;

 private:
  size_t n;
};

  class AnovaDimensionVarianceShareCutter : public AnovaCutter {
 public:
    AnovaDimensionVarianceShareCutter(double minCoveredVariance);

  AnovaResult cut(const SGridSample& input, const AnovaInfo& info) override;

 private:
  double minCoveredVariance;
  };

class AnovaComponentVarianceCutter : public AnovaCutter {
 public:
  AnovaComponentVarianceCutter(double minVariance);

  AnovaResult cut(const SGridSample& input, const AnovaInfo& info) override;

 private:
  double minVariance;
};

class AnovaReducer : public sgpp::base::Reducer<SGridSample, AnovaInfo, AnovaResult> {
 public:
  AnovaInfo evaluate(SGridSample& input) override;
};

}  // namespace base
}  // namespace sgpp

#endif
