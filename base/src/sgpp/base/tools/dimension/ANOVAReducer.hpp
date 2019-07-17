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

struct AnovaResult {
  std::vector<bool> activeDimensions;
  AnovaHelper::AnovaComponentVector activeComponents;
  double coveredVariance;
  size_t dimensions;

  TransformationFunction createTransformationFunction();
  GridSample<double> apply(GridSample<double>& sample);
};

class AnovaCutter : public sgpp::base::Cutter<GridSample<double>, AnovaInfo, AnovaResult> {};

class AnovaFixedCutter : public AnovaCutter {
 public:
  AnovaFixedCutter(size_t n);

  AnovaResult cut(const GridSample<double>& input, const AnovaInfo& info) override;

 private:
  size_t n;
};

  class AnovaDimensionVarianceShareCutter : public AnovaCutter {
 public:
    AnovaDimensionVarianceShareCutter(double minCoveredVariance);

  AnovaResult cut(const GridSample<double>& input, const AnovaInfo& info) override;

 private:
  double minCoveredVariance;
  };

class AnovaComponentVarianceCutter : public AnovaCutter {
 public:
  AnovaComponentVarianceCutter(double minVariance);

  AnovaResult cut(const GridSample<double>& input, const AnovaInfo& info) override;

 private:
  double minVariance;
};

class AnovaReducer : public sgpp::base::Reducer<GridSample<double>, AnovaInfo, AnovaResult> {
 public:
  AnovaInfo evaluate(GridSample<double>& input) override;
};

}  // namespace base
}  // namespace sgpp

#endif
