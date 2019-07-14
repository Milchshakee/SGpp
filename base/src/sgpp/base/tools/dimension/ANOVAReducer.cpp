#include <iostream>
#include "AnovaReducer.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "OperationAnova.hpp"
#include "DimReduction.hpp"

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
  size_t active = 0;
  for (size_t d = 0; d < activeDimensions.size(); d++) {
    if (activeDimensions[d]) {
      active++;
    }
  }

  bool hasZeroComponent = false;
  for (size_t i = 0; i < activeComponents.size(); i++) {
    if (activeComponents[i] == AnovaHelper::AnovaComponent(active, false)) {
      hasZeroComponent = true;
      break;
    }
  }

  AnovaHelper::AnovaComponentVector newComps(activeComponents);
  if (!hasZeroComponent) {
    newComps.emplace_back(AnovaHelper::AnovaComponent(active, false));
    }

  std::shared_ptr<Grid> newGrid(
      sgpp::base::Grid::createAnovaBoundaryGrid(active, activeComponents));
  newGrid->getGenerator().regular(const_cast<Grid&>(sample.getGrid()).getStorage().getMaxLevel());

    std::function<double(const DataVector&)> f = [sample](const DataVector& v) {
    return sample.getValue(v);
  };
  GridSample<double> newSample(newGrid, f);
  return std::move(newSample);
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
  OperationAnova op(const_cast<Grid&>(input.getGrid()).getStorage());
  sgpp::base::Sample<AnovaHelper::AnovaComponent, double> variances =
      op.calculateAnovaComponentVariances(v);

  return AnovaInfo{variances};
}