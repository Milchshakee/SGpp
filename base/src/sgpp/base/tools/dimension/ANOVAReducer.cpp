
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationEvalAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/tools/dimension/AnovaReducer.hpp>

sgpp::base::AnovaResult::AnovaResult(std::vector<bool>& ad, size_t d, const SGridSample& sample,
                                     double errorShare)
    : activeDimensions(ad), dimensions(d), errorShare(errorShare) {
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
  f = MatrixFunction(mat);

  // create reduced sample
  std::shared_ptr<Grid> newGrid(sgpp::base::Grid::createAnovaPrewaveletBoundaryGrid(active));
  newGrid->getGenerator().regular(const_cast<Grid&>(sample.getGrid()).getStorage().getMaxLevel());

  std::function<double(const DataVector&)> f = [this, sample](const DataVector& v) {
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
  reducedSample.hierarchise();
  eval = EvalFunction(reducedSample);
  originalFunction = EvalFunction(sample);
}

sgpp::base::ScalarFunction& sgpp::base::AnovaResult::getOriginalFunction() {
  return originalFunction;
}

sgpp::base::VectorFunction& sgpp::base::AnovaResult::getTransformationFunction() { return f; }

sgpp::base::ScalarFunction& sgpp::base::AnovaResult::getReducedFunction() { return eval; }

sgpp::base::SGridSample& sgpp::base::AnovaResult::getReducedOutput() { return reducedSample; }

sgpp::base::AnovaErrorRuleCutter::AnovaErrorRuleCutter(ErrorRule& r, double maxError)
    : ErrorRuleCutter<sgpp::base::SGridSample, sgpp::base::AnovaInfo, sgpp::base::AnovaResult>(
          r, maxError) {}

namespace {
void cutRec(double maxErrorShare, double* currentError, size_t* dim,
            sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double>& errors,
            std::vector<bool>& activeDimensions, size_t minDim = 1) {
  if (*dim <= minDim) {
    return;
  }

  std::vector<double> values(*dim);
  for (size_t d = 0; d < *dim; d++) {
    double er = 0;
    for (sgpp::base::AnovaBoundaryGrid::AnovaComponent& c : errors.getKeys()) {
      if (c[d]) {
        er += errors.getValue(c);
      }
    }
    values[d] = er;
  }

  auto min = std::min_element(values.begin(), values.end());
  size_t minIndex = min - values.begin();

  double additionalError = *min;
  if (((*currentError) + additionalError) < maxErrorShare) {
    // Update active components
    size_t removed = 0;
    for (size_t i = 0; i < errors.getSize(); i++) {
      if (errors.getKeys()[i][minIndex]) {
        errors.getKeys().erase(errors.getKeys().begin() + i - removed);
        removed++;
      }
    }

    // Remove dimension from remaining components
    for (sgpp::base::AnovaBoundaryGrid::AnovaComponent& c : errors.getKeys()) {
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

    (*currentError) += additionalError;
    (*dim) -= 1;
    cutRec(maxErrorShare, currentError, dim, errors, activeDimensions, minDim);
  }
}
}  // namespace

sgpp::base::AnovaFixedCutter::AnovaFixedCutter(size_t n)
    : FixedCutter<sgpp::base::SGridSample, sgpp::base::AnovaInfo, sgpp::base::AnovaResult>(n) {}

sgpp::base::AnovaResult sgpp::base::AnovaFixedCutter::cut(const SGridSample& input,
                                                          const AnovaInfo& info) {
  size_t dim = input.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double> errorShares =
      info.errorShares;
  double share = 1.0;
  cutRec(1.0, &share, &dim, errorShares, activeDimensions, n);

  return AnovaResult(activeDimensions, dim, input, share);
}

sgpp::base::AnovaResult sgpp::base::AnovaErrorRuleCutter::cut(const SGridSample& input,
                                                              const AnovaInfo& info) {
  size_t dim = input.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double> errorShares =
      info.errorShares;
  double share = 1.0;
  cutRec(maxError, &share, &dim, errorShares, activeDimensions);

  return AnovaResult(activeDimensions, dim, input, share);
}


sgpp::base::AnovaReducer::AnovaReducer(ErrorRule& rule) : rule(rule) {
}

sgpp::base::AnovaInfo sgpp::base::AnovaReducer::evaluate(SGridSample& input) {
  EvalFunction s(input);
  double sum = rule.calculateAbsoluteError(s);
  std::map<AnovaBoundaryGrid::AnovaComponent, double> errorShares;
  for (size_t i = 0; i < input.getSize(); i++) {
    GridPoint& gp = const_cast<Grid&>(input.getGrid()).getStorage().getPoint(i);
    AnovaBoundaryGrid::AnovaComponent c = AnovaBoundaryGrid::getAnovaComponentOfPoint(gp);
    if (errorShares.find(c) == errorShares.end()) {
      OperationEvalAnovaPrewaveletBoundary b(const_cast<Grid&>(input.getGrid()).getStorage(), c);
      EvalFunction f(input, b);
      double val = rule.calculateAbsoluteError(f) / sum;
      errorShares.emplace(c, val);
    }
  }
  return {errorShares};
}
