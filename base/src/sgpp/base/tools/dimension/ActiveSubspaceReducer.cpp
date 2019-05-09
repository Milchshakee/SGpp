#include "ActiveSubspaceReducer.hpp"
#include "Tools.hpp"

namespace sgpp {
namespace base {

ActiveSubspaceReducer::GivenGradient::GivenGradient(
    std::shared_ptr<sgpp::optimization::VectorFunction> gradient)
    : gradient(gradient) {}

sgpp::base::DataVector ActiveSubspaceReducer::GivenGradient::gradientAt(sgpp::base::DataVector& v) {
  sgpp::base::DataVector result(gradient->getNumberOfComponents());
  gradient->eval(v, result);
  return std::move(result);
}

ActiveSubspaceReducer::EigenValueCutoff::EigenValueCutoff(double minValue)
    : minEigenValue(minValue) {}

size_t ActiveSubspaceReducer::EigenValueCutoff::evaluate(const ActiveSubspaceInfo& info) {
  for (size_t d = 0; d < info.eigenValues.size(); ++d) {
    if (info.eigenValues[d] < minEigenValue) {
      return d + 1;
    }
  }
  return info.eigenValues.size();
}

ActiveSubspaceReducer::IntervalCutoff::IntervalCutoff(size_t bootstrapSamples) : bootstrapSamples(bootstrapSamples) {}

size_t ActiveSubspaceReducer::IntervalCutoff::evaluate(const ActiveSubspaceInfo& info) {
  size_t dimensions = info.eigenValues.size();
  std::mt19937_64 prng;
  std::uniform_int_distribution<size_t> dist(0, info.sampleMatrices.size() - 1);

  std::vector<std::pair<double, double>> eigenValueIntervals(dimensions);
  for (size_t d = 0; d < dimensions; ++d) {
    eigenValueIntervals[d].first = info.eigenValues[d];
    eigenValueIntervals[d].second = info.eigenValues[d];
  }

  for (size_t i = 0; i < bootstrapSamples; i++) {
    sgpp::base::DataMatrix bootstrapMatrix(info.eigenValues.size(), info.eigenValues.size());
    for (size_t j = 0; j < info.sampleMatrices.size(); i++) {
      size_t l = dist(prng);
      bootstrapMatrix.add(info.sampleMatrices[l]);
    }
    bootstrapMatrix.mult(1.0 / static_cast<double>(bootstrapSamples));

    sgpp::base::DataMatrix bootstrapEigenVectorMatrix(dimensions, dimensions);
    sgpp::base::DataVector bootstrapEigenValues(dimensions);
    Tools::svd(bootstrapMatrix, bootstrapEigenVectorMatrix, bootstrapEigenValues);
    for (size_t d = 0; d < dimensions; ++d) {
      double e = bootstrapEigenValues[d];
      if (e < eigenValueIntervals[d].first) {
        eigenValueIntervals[d].first = e;
      }
      if (e > eigenValueIntervals[d].second) {
        eigenValueIntervals[d].second = e;
      }
    }
  }

  double max = 0;
  size_t cutoff = 0;
  for (size_t d = 0; d < dimensions; ++d) {
    double size = eigenValueIntervals[d].second - eigenValueIntervals[d].first;
    if (size > max) {
      max = size;
      cutoff = d + 1;
    }
  }
  return cutoff;
}

ActiveSubspaceReducer::ActiveSubspaceReducer(size_t samples, std::shared_ptr<Gradient> gradient,
  std::shared_ptr<VectorDistribution> distribution,
    std::shared_ptr<CutoffCriterion<ActiveSubspaceInfo>> cutoff) : FunctionReducer<sgpp::base::ActiveSubspaceInfo>(cutoff), samples(samples), gradient(gradient), distribution(distribution) {}


void ActiveSubspaceReducer::evaluateFunction(
    sgpp::optimization::ScalarFunction& input, ActiveSubspaceInfo& out) {
  size_t dimensions = input.getNumberOfParameters();
  sgpp::base::DataMatrix matrix(dimensions, dimensions);
  out.sampleMatrices = std::vector<DataMatrix>(samples, sgpp::base::DataMatrix(dimensions, dimensions));
  for (size_t i = 0; i < samples; ++i) {
    sgpp::base::DataVector sample = distribution->operator()();
    sgpp::base::DataVector sampleGradient = gradient->gradientAt(sample);
    for (size_t d = 0; d < dimensions; ++d) {
      sgpp::base::DataVector col = sampleGradient;
      col.mult(sampleGradient[d]);
      out.sampleMatrices[i].setColumn(d, col);
    }
    matrix.add(out.sampleMatrices[i]);
  }
  matrix.mult(1.0 / static_cast<double>(samples));

  out.eigenVectors = sgpp::base::DataMatrix(dimensions, dimensions);
  out.eigenValues = sgpp::base::DataVector(dimensions);
  Tools::svd(matrix, out.eigenVectors, out.eigenValues);
}

std::unique_ptr<sgpp::optimization::ScalarFunction> ActiveSubspaceReducer::reduce(
    sgpp::optimization::ScalarFunction& input, size_t n, const ActiveSubspaceInfo& info) {
  DataMatrix m = info.eigenVectors;
  m.resizeRowsCols(input.getNumberOfParameters(), n);
  std::unique_ptr<optimization::ScalarFunction> ptr;
  input.clone(ptr);
  return std::make_unique<ReducedFunction>(std::move(ptr), m);
}

}  // namespace base
}  // namespace sgpp
