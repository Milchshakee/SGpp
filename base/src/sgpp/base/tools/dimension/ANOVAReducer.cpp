
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/dimension/AnovaReducer.hpp>
#include <sgpp/base/tools/dimension/OperationAnova.hpp>

sgpp::base::AnovaResult::AnovaResult(std::vector<bool>& ad, size_t d, const SGridSample& sample)
    : activeDimensions(ad), dimensions(d) {
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
  std::shared_ptr<Grid> newGrid(sgpp::base::Grid::createAnovaBoundaryGrid(active));
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

double calcError(sgpp::base::AnovaResult& r, const sgpp::base::SGridSample& sample,
                 sgpp::base::ErrorRule& rule, size_t dims,
                 size_t removedDim, std::vector<bool>& ad) {
  std::vector<bool> activeCopy = ad;
  ad[removedDim] = false;
  r = sgpp::base::AnovaResult(activeCopy, dims - 1, sample);
  return r.calculateRelativeError(rule);
}

sgpp::base::AnovaResult cutRec(const sgpp::base::SGridSample& sample, double maxError,
                               sgpp::base::ErrorRule& rule, size_t originalDim, size_t* dim,
                               std::vector<bool>& activeDimensions, size_t minDim = 1) {
  sgpp::base::AnovaResult r;
  double minError = std::numeric_limits<double>::infinity();
  size_t min = 0;
  for (size_t d = 0; d < originalDim; d++) {
    if (activeDimensions[d]) {
      double error = calcError(r, sample, rule, *dim, d, activeDimensions);
      if (error < minError) {
        minError = error;
        min = d;
      }
    }
  }

  if (minError <= maxError) {
    if (*dim <= minDim) {
      return r;
    }

    activeDimensions[min] = false;
    (*dim) -= 1;
    cutRec(sample, maxError, rule, originalDim, dim, activeDimensions, minDim);
  }
}
}  // namespace

sgpp::base::AnovaFixedCutter::AnovaFixedCutter(ErrorRule& r, size_t n)
    : FixedCutter<sgpp::base::SGridSample, sgpp::base::AnovaInfo, sgpp::base::AnovaResult>(r, n) {}

sgpp::base::AnovaResult sgpp::base::AnovaFixedCutter::cut(const SGridSample& input,
                                                          const AnovaInfo& info) {
  size_t dim = input.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  cutRec(input, 1.0, r, input.getDimensions(), &dim, activeDimensions, n);

  return AnovaResult(activeDimensions, dim, input);
}

sgpp::base::AnovaResult sgpp::base::AnovaErrorRuleCutter::cut(const SGridSample& input,
                                                              const AnovaInfo& info) {
  size_t dim = input.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  cutRec(input, maxError, r, input.getDimensions(), &dim, activeDimensions);

  return AnovaResult(activeDimensions, dim, input);
}

sgpp::base::AnovaInfo sgpp::base::AnovaReducer::evaluate(SGridSample& input) { return AnovaInfo{}; }