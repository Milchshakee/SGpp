#include "AnovaReducer.hpp"
#include "DimReduction.hpp"
#include "OperationAnova.hpp"
#include "sgpp/base/grid/Grid.hpp"

sgpp::base::AnovaComponentVarianceCutter::AnovaComponentVarianceCutter(double minVariance)
    : minVariance(minVariance) {}

sgpp::base::AnovaResult::AnovaResult(std::vector<bool>& activeDimensions,
                                     AnovaHelper::AnovaComponentVector& activeComponents,
                                     double coveredVariance, size_t dimensions, SGridSample sample)
    : activeDimensions(activeDimensions),
      activeComponents(activeComponents),
      coveredVariance(coveredVariance),
      dimensions(dimensions) {
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
  reducedSample = SGridSample(newGrid, f);
}

double sgpp::base::AnovaResult::calcMcL2Error(optimization::ScalarFunction& func, size_t paths,
                                              uint64_t seed) {
  std::mt19937_64 rand(seed);
  std::uniform_real_distribution<double> dist(0, 1);
  size_t dim = func.getNumberOfParameters();
  BoundingBox& boundingBox = const_cast<Grid&>(reducedSample.getGrid()).getBoundingBox();

  sgpp::base::DataVector point(dim);
  double res = 0;

  for (size_t i = 0; i < paths; i++) {
    for (size_t d = 0; d < dim; d++) {
      point[d] = boundingBox.transformPointToBoundingBox(d, dist(rand));
    }
    double val = func.eval(point);
    res += pow(val - reducedSample.eval(point), 2);
  }

  double determinant = 1.0;
  for (size_t d = 0; d < dim; d++) {
    determinant *= boundingBox.getIntervalWidth(d);
  }

  return sqrt(res / static_cast<double>(paths) * determinant);
}

sgpp::base::TransformationFunction& sgpp::base::AnovaResult::getTransformationFunction() {
  return f;
}

sgpp::base::SGridSample& sgpp::base::AnovaResult::getReducedOutput() { return reducedSample; }

sgpp::base::AnovaResult sgpp::base::AnovaComponentVarianceCutter::cut(const SGridSample& input,
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
  for (AnovaHelper::AnovaComponent c : comps) {
    coveredVar += info.variances.getValue(c);
  }

  return AnovaResult(activeDimensions, comps, coveredVar, dim, input);
}

sgpp::base::AnovaDimensionVarianceShareCutter::AnovaDimensionVarianceShareCutter(
    double minCoveredVariance)
    : minCoveredVariance(minCoveredVariance) {}

void cutRec(double minShare, double* currentVar, double totalVar, size_t* dim,
            const sgpp::base::AnovaInfo& info, sgpp::base::AnovaHelper::AnovaComponentVector& comps,
            std::vector<bool>& activeDimensions, size_t minDim = 1) {
  if (*dim < minDim) {
    return;
    }

  std::vector<double> vars(*dim);
  for (size_t d = 0; d < *dim; d++) {
    double var = 0;
    for (sgpp::base::AnovaHelper::AnovaComponent c : comps) {
      if (c[d]) {
        var += info.variances.getValue(c);
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
    for (sgpp::base::AnovaHelper::AnovaComponent& c : comps) {
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
    cutRec(minShare, currentVar, totalVar, dim, info, comps, activeDimensions, minDim);
  }
}

sgpp::base::AnovaFixedCutter::AnovaFixedCutter(size_t n) : n(n) {}

sgpp::base::AnovaResult sgpp::base::AnovaFixedCutter::cut(const SGridSample& input,
                                                          const AnovaInfo& info) {
  AnovaHelper::AnovaComponentVector comps(info.variances.getKeys());
  double var = info.totalVariance;
  size_t dim = input.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  cutRec(0.0, &var, info.totalVariance, &dim, info, comps, activeDimensions, n);

  return AnovaResult(activeDimensions, comps, var, dim, input);
}

sgpp::base::AnovaResult sgpp::base::AnovaDimensionVarianceShareCutter::cut(const SGridSample& input,
                                                                           const AnovaInfo& info) {
  AnovaHelper::AnovaComponentVector comps(info.variances.getKeys());
  double var = info.totalVariance;
  size_t dim = input.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  cutRec(minCoveredVariance, &var, info.totalVariance, &dim, info, comps, activeDimensions);

  return AnovaResult(activeDimensions, comps, var, dim, input);
}

sgpp::base::AnovaInfo sgpp::base::AnovaReducer::evaluate(SGridSample& input) {
  DataVector v(input.getValues());
  OperationAnova op(const_cast<Grid&>(input.getGrid()).getStorage());
  sgpp::base::Sample<AnovaHelper::AnovaComponent, double> variances =
      op.calculateAnovaComponentVariances(v);
  double sum = std::accumulate(variances.getValues().begin(), variances.getValues().end(), 0);
  return AnovaInfo{variances, sum};
}