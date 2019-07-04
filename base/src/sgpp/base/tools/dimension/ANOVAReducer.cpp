#include <iostream>
#include "AnovaReducer.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"

sgpp::base::AnovaFixedCutter::AnovaFixedCutter(size_t n) : n(n) {}

sgpp::base::AnovaResult sgpp::base::AnovaFixedCutter::cut(const GridSample<double>& input,
                                                          const AnovaInfo& info) {}

sgpp::base::TransformationFunction sgpp::base::AnovaResult::createTransformationFunction() {
  DataMatrix mat(0, activeDimensions.size());
  for (size_t d = 0; d < activeDimensions.size(); d++) {
    if (activeDimensions[d]) {
      DataVector v(activeDimensions.size(), 0);
      v[d] = 1;
      mat.appendRow(v);
    }
  }
  return TransformationFunction(mat);
}

sgpp::base::GridSample<double> sgpp::base::AnovaResult::apply(GridSample<double>& sample) {
  if (activeDimensions.empty()) {
    //TODO
    }


  size_t active = 0;
  for (size_t d = 0; d < activeDimensions.size(); d++) {
    if (activeDimensions[d]) {
      active++;
    }
  }

  std::shared_ptr<Grid> newGrid(
      sgpp::base::Grid::createAnovaBoundaryGrid(active, activeComponents));
  newGrid->getGenerator().regular(const_cast<Grid&>(sample.getGrid()).getStorage().getMaxLevel());

  TransformationFunction t = createTransformationFunction();
  std::map<DataVector, std::tuple<size_t, double>> collapsingPoints;
  GridDistribution dist(*newGrid);
  for (const DataVector& v : dist.getVectors()) {
    collapsingPoints.emplace(v, std::tuple<size_t, double>{0, 0.0});
  }

  for (const DataVector& v : sample.getKeys()) {
    DataVector newV;
    t.eval(v, newV);
    std::tuple<size_t, double>& t = collapsingPoints.at(newV);
    std::get<0>(t)++;
    std::get<1>(t) += sample.getValue(v);
  }

  std::function<double(const DataVector&)> f = [collapsingPoints](const DataVector& v) {
    auto t = collapsingPoints.at(v);
    return std::get<1>(t) / static_cast<double>(std::get<0>(t));
  };
  GridSample<double> s(newGrid, f);
  return std::move(s);
}

sgpp::base::AnovaVarianceCutter::AnovaVarianceCutter(double minVariance)
    : minVariance(minVariance) {}

sgpp::base::AnovaResult sgpp::base::AnovaVarianceCutter::cut(const GridSample<double>& input,
                                                             const AnovaInfo& info) {
  AnovaHelper::AnovaComponentVector comps(info.variances.getKeys());
  std::vector<bool> activeDimensions(input.getDimensions(), false);
  size_t removed = 0;
  for (size_t i = 0; i < info.variances.getSize(); i++) {
    if (info.variances.getValues()[i] < minVariance) {
      comps.erase(comps.begin() + i - removed);
      removed++;
    } else {
      for (size_t d = 0; d < input.getGrid().getDimension(); d++) {
        if (info.variances.getKeys()[i][d]) {
          activeDimensions[d] = true;
        }
      }
    }
  }

  for (size_t d = 0; d < input.getGrid().getDimension(); d++) {
    if (!activeDimensions[d]) {
      for (size_t i = 0; i < info.variances.getSize(); i++) {
        comps[i].erase(comps[i].begin() + d);
      }
    }
  }

  return AnovaResult{activeDimensions, comps};
}

sgpp::base::AnovaInfo sgpp::base::AnovaReducer::evaluate(GridSample<double>& input) {
  DataVector v(input.getValues());
  std::unique_ptr<sgpp::base::OperationHierarchisation>(
      sgpp::op_factory::createOperationHierarchisation(const_cast<Grid&>(input.getGrid())))
      ->doHierarchisation(v);
  std::unique_ptr<sgpp::base::OperationAnova> op(
      sgpp::op_factory::createOperationAnova(const_cast<Grid&>(input.getGrid())));
  sgpp::base::Sample<AnovaHelper::AnovaComponent, double> variances =
      op->calculateAnovaOrderVariances(v);

  return AnovaInfo{variances};
}