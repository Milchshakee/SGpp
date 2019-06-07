// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ANOVAREDUCER_HPP
#define ANOVAREDUCER_HPP

#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include "sgpp/base/tools/Sample.hpp"
#include "sgpp/base/grid/type/AnovaBoundaryGrid.hpp"

namespace sgpp {
namespace base {

struct AnovaInfo
{
  sgpp::base::Sample<sgpp::base::AnovaComponent, double> variances;
};

struct AnovaResult
{
  AnovaComponentVector activeComponents;

  std::unique_ptr<sgpp::base::AnovaBoundaryGrid> apply(sgpp::base::AnovaBoundaryGrid& grid);
};

    class AnovaCutter : public sgpp::base::Cutter<AnovaInfo, GridSample<double>, AnovaResult> {};

  
  class AnovaVarianceCutter : public sgpp::base::Cutter<AnovaInfo, GridSample<double>, AnovaResult> {
 public:
    AnovaVarianceCutter(double maxVariance);

  AnovaResult cut(const AnovaInfo& input, const GridSample<double>& info) override;

 private:
  double maxVariance;
  };

class AnovaReducer : public sgpp::base::Reducer<GridSample<double>,AnovaInfo,AnovaResult> {
 public:

  AnovaReducer();

  void evaluate(GridSample<double>& input, AnovaInfo& out) override;
};

}  // namespace base
}  // namespace sgpp

#endif
