#include <iostream>
#include "AnovaReducer.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"


sgpp::base::TransformationFunction sgpp::base::AnovaResult::createTransformationFunction() {
  DataMatrix mat = DataMatrix();
  for (size_t d = 0; d < activeDimensions.size(); d++) {
    if (activeDimensions[d]) {
      DataVector v(activeDimensions.size(), 0);
      v[d] = 1;
      mat.appendRow(v);
    }
  }
  return TransformationFunction(mat);
}

sgpp::base::GridSample<double> sgpp::base::AnovaResult::createReducedGridSample(
  GridSample<double>& sample) {
  size_t active = 0;
  for (size_t d = 0; d < activeDimensions.size(); d++) {
    if (activeDimensions[d]) {
      active++;
    }
  }

  std::shared_ptr<Grid> newGrid(
      sgpp::base::Grid::createAnovaBoundaryGrid(active, activeComponents));
  TransformationFunction t = createTransformationFunction();
  std::map<DataVector,std::tuple<size_t, double>> collapsingPoints;
  for (const DataVector& v : GridDistribution(*newGrid).getVectors()) {
    collapsingPoints.emplace(v, std::tuple<size_t,double> {0, 0.0});
  }

  for (const DataVector& v : sample.getKeys()) {
    DataVector newV;
    t.eval(v, newV);
    std::tuple<size_t, double>& t = collapsingPoints.at(newV);
    std::get<0>(t)++;
    std::get<1>(t)+= sample.getValue(v);
  }

  std::vector<double> newValues(newGrid->getSize());
  std::function<double(const DataVector&)> f = [collapsingPoints](const DataVector& v){
    auto t = collapsingPoints.at(v);
    return static_cast<double>(std::get<1>(t)) / std::get<0>(t);
  };
  GridSample<double> s(newGrid, f);
  return std::move(s);
}

sgpp::base::AnovaVarianceCutter::AnovaVarianceCutter(double maxVariance) : maxVariance(maxVariance) {
}

sgpp::base::AnovaResult sgpp::base::AnovaVarianceCutter::cut(
    const GridSample<double>& input, const AnovaInfo& info)
{
  AnovaComponentVector comps(info.variances.getKeys());
  size_t removed = 0;
  std::vector<bool> activeDimensions(input.getDimensions(), true);

  out:
  for (size_t d = 0; d < input.getDimensions(); d++) {
    for (const AnovaComponent& comp : info.variances.getKeys()) {
      if (comp[d]) {
        goto out;
        }
      }

        for (AnovaComponent comp : comps) {
        comps.erase(comps.begin() + d - removed);
      }
        removed++;
      activeDimensions[d] = false;
  }
  return AnovaResult{activeDimensions, comps};
  }

sgpp::base::AnovaInfo sgpp::base::AnovaReducer::evaluate(GridSample<double>& input) {
    DataVector v(input.getValues());
    std::unique_ptr<sgpp::base::OperationHierarchisation>(
        sgpp::op_factory::createOperationHierarchisation(const_cast<Grid&>(input.getGrid())))
        ->doHierarchisation(v);
  std::unique_ptr<sgpp::base::OperationAnova> op(sgpp::op_factory::createOperationAnova(const_cast<Grid&>(input.getGrid())));  
    sgpp::base::Sample<sgpp::base::AnovaComponent, double> variances =
      op->calculateAnovaOrderVariances(v);

  return AnovaInfo{variances};
}