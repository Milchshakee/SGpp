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
  sgpp::base::Sample<AnovaBoundaryGrid::AnovaComponent, double> variances;
  double totalVariance;
};

class AnovaResult : public Result<SGridSample> {
public:
  std::vector<bool> activeDimensions;
 size_t dimensions;

  AnovaResult(std::vector<bool>& activeDimensions, double coveredVariance,
             size_t dimensions, SGridSample sample);

  VectorFunction& getTransformationFunction();
  SGridSample& getReducedOutput();
  ScalarFunction& getReducedFunction() override;
  double getCoveredVariance();

 private:
  double coveredVariance;
  EvalFunction eval;
  MatrixFunction f;
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

  class AnovaVarianceShareCutter : public AnovaCutter {
 public:
    AnovaVarianceShareCutter(double minCoveredVariance);

  AnovaResult cut(const SGridSample& input, const AnovaInfo& info) override;

 private:
  double minCoveredVariance;
  };

class AnovaReducer : public sgpp::base::Reducer<SGridSample, AnovaInfo, AnovaResult> {
 public:
  AnovaInfo evaluate(SGridSample& input) override;
};

}  // namespace base
}  // namespace sgpp

#endif
