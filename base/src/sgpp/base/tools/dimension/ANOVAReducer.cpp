#include "AnovaReducer.hpp"
#include "DimReduction.hpp"
#include "OperationAnova.hpp"
#include "sgpp/base/grid/Grid.hpp"

sgpp::base::AnovaComponentVarianceCutter::AnovaComponentVarianceCutter(double minVariance)
    : minVariance(minVariance) {}

sgpp::base::AnovaResult::AnovaResult(std::vector<bool>& ad,
                                     double cv, size_t d, SGridSample sample)
    : activeDimensions(ad),
      coveredVariance(cv),
      dimensions(d) {
  size_t active = 0;
  for (size_t d = 0; d < activeDimensions.size(); d++) {
    if (activeDimensions[d]) {
      active++;
    }
  }
  if (active == 0) {
    activeDimensions[0] = true;
    active++;
    dimensions++;
  }

  // create transformation function

  DataMatrix mat(0, activeDimensions.size());
  for (size_t d = 0; d < activeDimensions.size(); d++) {
    if (activeDimensions[d]) {
      DataVector v(activeDimensions.size(), 0);
      v[d] = 1;
      mat.appendRow(v);
    }
  }
  f = TransformationFunction(mat);

  // create reduced sample
  std::shared_ptr<Grid> newGrid(
      sgpp::base::Grid::createAnovaBoundaryGrid(active));
  newGrid->getGenerator().regular(const_cast<Grid&>(sample.getGrid()).getStorage().getMaxLevel());

  std::function<double(const DataVector&)> f = [this,sample](const DataVector& v) {
    DataVector newV(sample.getDimensions());
    size_t counter = 0;
    for (size_t d = 0; d < activeDimensions.size(); d++) {
      if (activeDimensions[d]) {
        newV[d] = v[counter];
        counter++;
      } else {
        newV[d] = 0;
        }
    }
    return sample.getValue(newV);
  };
  reducedSample = SGridSample(newGrid, f);
  eval = EvalFunction(reducedSample);
}

sgpp::base::TransformationFunction& sgpp::base::AnovaResult::getTransformationFunction() {
  return f;
}

sgpp::optimization::ScalarFunction& sgpp::base::AnovaResult::getReducedFunction() { return eval; }

sgpp::base::SGridSample& sgpp::base::AnovaResult::getReducedOutput() { return reducedSample; }

sgpp::base::AnovaResult sgpp::base::AnovaComponentVarianceCutter::cut(const SGridSample& input,
                                                                      const AnovaInfo& info) {
  AnovaBoundaryGrid::AnovaComponentVector comps(info.variances.getKeys());
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

  size_t dim = 0;
  for (size_t d = 0; d < input.getGrid().getDimension(); d++) {
    if (!activeDimensions[d]) {
      for (size_t i = 0; i < info.variances.getSize(); i++) {
        comps[i].erase(comps[i].begin() + d);
      }
    } else {
      dim++;
    }
  }

  double coveredVar = 0;
  for (AnovaBoundaryGrid::AnovaComponent c : comps) {
    coveredVar += info.variances.getValue(c);
  }

  return AnovaResult(activeDimensions, coveredVar, dim, input);
}

sgpp::base::AnovaDimensionVarianceShareCutter::AnovaDimensionVarianceShareCutter(
    double minCoveredVariance)
    : minCoveredVariance(minCoveredVariance) {}

void cutRec(double minShare, double* currentVar, double totalVar, size_t* dim,
            sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double>& variances,
            std::vector<bool>& activeDimensions, size_t minDim = 1) {
  if (*dim <= minDim) {
    return;
  }

  std::vector<sgpp::base::AnovaBoundaryGrid::AnovaComponent>& comps =
      const_cast<std::vector<sgpp::base::AnovaBoundaryGrid::AnovaComponent>&>(variances.getKeys());

  std::vector<double> vars(*dim);
  for (size_t d = 0; d < *dim; d++) {
    double var = 0;
    for (sgpp::base::AnovaBoundaryGrid::AnovaComponent& c : comps) {
      if (c[d]) {
        var += variances.getValue(c);
      }
    }
    vars[d] = var;
  }

  auto min = std::min_element(vars.begin(), vars.end());
  size_t minIndex = min - vars.begin();

  double removedVar = *min;
  if ((((*currentVar) - removedVar) / totalVar) >= minShare) {
    // Update active components
    size_t removed = 0;
    for (size_t i = 0; i < comps.size(); i++) {
      if (comps[i][minIndex]) {
        comps.erase(comps.begin() + i - removed);
        removed++;
      }
    }

    // Remove dimension from remaining components
    for (sgpp::base::AnovaBoundaryGrid::AnovaComponent& c : comps) {
      c.erase(c.begin() + minIndex);
    }

    // Update active dimensions
    size_t counter = 0;
    for (size_t d = 0; d < activeDimensions.size(); d++) {
      if (activeDimensions[d]) {
        if (counter == minIndex) {
          activeDimensions[d] = false;
        }
        counter++;
      }
    }

    (*currentVar) -= removedVar;
    (*dim) -= 1;
    cutRec(minShare, currentVar, totalVar, dim, variances, activeDimensions, minDim);
  }
}

sgpp::base::AnovaFixedCutter::AnovaFixedCutter(size_t n) : n(n) {}

sgpp::base::AnovaResult sgpp::base::AnovaFixedCutter::cut(const SGridSample& input,
                                                          const AnovaInfo& info) {
  double var = info.totalVariance;
  size_t dim = input.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double> variances = info.variances;
  cutRec(0.0, &var, info.totalVariance, &dim, variances, activeDimensions, n);

  return AnovaResult(
      activeDimensions, var,
      dim, input);
}

sgpp::base::AnovaResult sgpp::base::AnovaDimensionVarianceShareCutter::cut(const SGridSample& input,
                                                                           const AnovaInfo& info) {
  double var = info.totalVariance;
  size_t dim = input.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double> variances = info.variances;
  cutRec(minCoveredVariance, &var, info.totalVariance, &dim, variances, activeDimensions);

  return AnovaResult(
      activeDimensions,
      var / info.totalVariance, dim, input);
}

sgpp::base::AnovaInfo sgpp::base::AnovaReducer::evaluate(SGridSample& input) {
  DataVector v(input.getValues());
  OperationAnova op(const_cast<Grid&>(input.getGrid()).getStorage());
  sgpp::base::Sample<AnovaBoundaryGrid::AnovaComponent, double> variances =
      op.calculateAnovaComponentVariances(v);
  double sum = std::accumulate(variances.getValues().begin(), variances.getValues().end(), 0);
  return AnovaInfo{variances, sum};
}