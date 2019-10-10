
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationEvalAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/tools/dimension/AnovaReducer.hpp>
#include <sgpp/base/operation/hash/OperationEvalAnovaLinearBoundary.hpp>

sgpp::base::AnovaResult::AnovaResult(std::vector<bool>& ad, size_t d, const AnovaInput& input,
                                     double converedVariance)
    : activeDimensions(ad), dimensions(d), coveredVariance(converedVariance), originalFunction(&input.originalFunction) {
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
  transformation = MatrixFunction(mat);
  DataMatrix inverse = mat;
  inverse.transpose();

  // create reduced sample
  std::shared_ptr<Grid> newGrid(sgpp::base::Grid::createLinearBoundaryGrid(active));
  size_t l = const_cast<Grid&>(input.sample.getGrid()).getStorage().getMaxLevel();
  newGrid->getGenerator().regular(l);

  EvalFunction sampleEval(input.sample);
  std::function<double(const DataVector&)> f = [this, &sampleEval, &inverse, &input](const DataVector& v) {
    DataVector newV(input.sample.getDimensions());
    inverse.mult(v, newV);
    return sampleEval.eval(newV);
  };
  reducedSample = SGridSample(newGrid, f);
  reducedSample.hierarchise();
  eval = EvalFunction(reducedSample);
}

sgpp::base::ScalarFunction& sgpp::base::AnovaResult::getOriginalFunction() {
  return *originalFunction;
}

sgpp::base::VectorFunction& sgpp::base::AnovaResult::getTransformationFunction() { return transformation; }

sgpp::base::ScalarFunction& sgpp::base::AnovaResult::getReducedFunctionSurrogate() { return eval; }

sgpp::base::SGridSample& sgpp::base::AnovaResult::getReducedOutput() { return reducedSample; }

sgpp::base::AnovaErrorRuleCutter::AnovaErrorRuleCutter(ErrorRule& r, double maxError)
    : ErrorRuleCutter<sgpp::base::AnovaInput, sgpp::base::AnovaInfo, sgpp::base::AnovaResult>(
          r, maxError) {}

namespace {
void cutRec(double minVarianceShare, double* currentVarianceShare, size_t* dim,
            sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double>& variances,
            std::vector<bool>& activeDimensions, size_t minDim = 1) {
  if (*dim <= minDim) {
    return;
  }

  std::vector<double> values(*dim);
  for (size_t d = 0; d < *dim; d++) {
    double er = 0;
    for (sgpp::base::AnovaBoundaryGrid::AnovaComponent& c : variances.getKeys()) {
      if (c[d]) {
        er += variances.getValue(c);
      }
    }
    values[d] = er;
  }

  auto min = std::min_element(values.begin(), values.end());
  size_t minIndex = min - values.begin();

  double removedVariance = *min;
  if (((*currentVarianceShare) - removedVariance) >= minVarianceShare) {
    // Update active components
    size_t removed = 0;
    for (size_t i = 0; i < variances.getSize(); i++) {
      if (variances.getKeys()[i - removed][minIndex]) {
        variances.getKeys().erase(variances.getKeys().begin() + i - removed);
        removed++;
      }
    }

    // Remove dimension from remaining components
    for (sgpp::base::AnovaBoundaryGrid::AnovaComponent& c : variances.getKeys()) {
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

    (*currentVarianceShare) -= removedVariance;
    (*dim) -= 1;
    cutRec(minVarianceShare, currentVarianceShare, dim, variances, activeDimensions, minDim);
  }
}
}  // namespace

sgpp::base::AnovaFixedCutter::AnovaFixedCutter(size_t n)
    : FixedCutter<sgpp::base::AnovaInput, sgpp::base::AnovaInfo, sgpp::base::AnovaResult>(n) {}

sgpp::base::AnovaResult sgpp::base::AnovaFixedCutter::cut(const AnovaInput& input,
                                                          const AnovaInfo& info) {
  size_t dim = input.sample.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double> varianceShares =
      info.variances;
  double share = 1.0;
  cutRec(0.0, &share, &dim, varianceShares, activeDimensions, n);

  return AnovaResult(activeDimensions, dim, input, share);
}

sgpp::base::AnovaResult sgpp::base::AnovaErrorRuleCutter::cut(const AnovaInput& input,
                                                              const AnovaInfo& info) {
  size_t dim = input.sample.getDimensions();
  std::vector<bool> activeDimensions(dim, true);
  sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double> errorShares =
      info.variances;
  double share = 1.0;
  cutRec(maxError, &share, &dim, errorShares, activeDimensions);

  return AnovaResult(activeDimensions, dim, input, share);
}


sgpp::base::AnovaReducer::AnovaReducer(ErrorRule& rule) : rule(rule) {
}

double evalComponent(sgpp::base::ErrorRule& rule, sgpp::base::AnovaInput& input,
                     sgpp::base::AnovaBoundaryGrid::AnovaComponent& comp)
{
  sgpp::base::OperationEvalAnovaLinearBoundary lin(const_cast<sgpp::base::Grid&>(input.sample.getGrid()).getStorage(), comp);
  sgpp::base::OperationEvalAnovaPrewaveletBoundary pre(const_cast<sgpp::base::Grid&>(input.sample.getGrid()).getStorage(),
                                           comp);
  sgpp::base::EvalFunction f;
  if (const_cast<sgpp::base::Grid&>(input.sample.getGrid()).getType() == sgpp::base::GridType::AnovaLinearBoundary) {
    f = sgpp::base::EvalFunction(input.sample, lin);
  } else {
    f = sgpp::base::EvalFunction(input.sample, pre);
  }
  double val = rule.calculateAbsoluteError(f);
  return val;
}


void iterateComponentsRec(
    std::map<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double>& errorShares, sgpp::base::ErrorRule& rule,
    sgpp::base::AnovaInput& input,
                                             sgpp::base::AnovaBoundaryGrid::AnovaComponent& comp, size_t dim, size_t maxDim) {
  if (dim == maxDim) {
    double v = evalComponent(rule, input, comp);
    errorShares.emplace(comp, v);
    return;
  }

  sgpp::base::AnovaBoundaryGrid::AnovaComponent copy = comp;
  copy[dim] = false;
  iterateComponentsRec(errorShares, rule, input, copy, dim + 1, maxDim);
  copy[dim] = true;
  iterateComponentsRec(errorShares, rule, input, copy, dim + 1, maxDim);
}

sgpp::base::AnovaInfo sgpp::base::AnovaReducer::evaluate(AnovaInput& input) {
  EvalFunction s(input.sample);
  AnovaBoundaryGrid::AnovaComponent start(input.sample.getDimensions(), false);
  std::map<AnovaBoundaryGrid::AnovaComponent, double> errorShares;
  iterateComponentsRec(
      errorShares, rule, input, start, 0,
      input.sample.getDimensions());

  double sum = 0;
  for (auto& x : errorShares) {
    sum += x.second;
    }

  for (auto& x : errorShares) {
      x.second /= sum;
  }


  return {errorShares};
}
