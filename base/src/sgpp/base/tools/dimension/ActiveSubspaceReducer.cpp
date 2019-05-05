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

ActiveSubspaceReducer::FixedCutoff::FixedCutoff(size_t n) : n(n) {}

size_t ActiveSubspaceReducer::FixedCutoff::evaluate(
    const sgpp::base::DataMatrix& eigenVectors, const sgpp::base::DataVector& eigenValues,
    const std::vector<sgpp::base::DataMatrix>& sampleMatrices) {
  return n;
}

size_t ActiveSubspaceReducer::EigenValueCutoff::evaluate(
    const sgpp::base::DataMatrix& eigenVectors, const sgpp::base::DataVector& eigenValues,
    const std::vector<sgpp::base::DataMatrix>& sampleMatrices) {
  for (size_t d = 0; d < eigenValues.size(); ++d) {
    if (eigenValues[d] < minEigenValue) {
      return d + 1;
    }
  }
  return eigenValues.size();
}

ActiveSubspaceReducer::IntervalCutoff::IntervalCutoff(size_t bootstrapSamples) : bootstrapSamples(bootstrapSamples) {}

size_t ActiveSubspaceReducer::IntervalCutoff::evaluate(
    const sgpp::base::DataMatrix& eigenVectors, const sgpp::base::DataVector& eigenValues,
    const std::vector<sgpp::base::DataMatrix>& sampleMatrices) {
  size_t dimensions = eigenValues.size();
  std::mt19937_64 prng;
  std::uniform_int_distribution<size_t> dist(0, sampleMatrices.size() - 1);

  std::vector<std::pair<double, double>> eigenValueIntervals(dimensions);
  for (size_t d = 0; d < dimensions; ++d) {
    eigenValueIntervals[d].first = eigenValues[d];
    eigenValueIntervals[d].second = eigenValues[d];
  }

  for (size_t i = 0; i < bootstrapSamples; i++) {
    sgpp::base::DataMatrix bootstrapMatrix(eigenValues.size(), eigenValues.size());
    for (size_t j = 0; j < sampleMatrices.size(); i++) {
      size_t l = dist(prng);
      bootstrapMatrix.add(sampleMatrices[l]);
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

ActiveSubspaceReducer::EigenValueCutoff::EigenValueCutoff(double minValue)
    : minEigenValue(minValue) {}

ActiveSubspaceReducer::ActiveSubspaceReducer(size_t samples, std::shared_ptr<Gradient> gradient,
                                             std::shared_ptr<VectorDistribution> distribution,
                                             std::shared_ptr<CutoffCriterion> cutoff)
    : samples(samples), gradient(gradient), distribution(distribution), cutoff(cutoff) {}

std::unique_ptr<sgpp::optimization::ScalarFunction> ActiveSubspaceReducer::reduceFunction(
    sgpp::optimization::ScalarFunction& input) {
  size_t dimensions = input.getNumberOfParameters();
  sgpp::base::DataMatrix matrix(dimensions, dimensions);
  std::vector<sgpp::base::DataMatrix> sampleMatrices(
      samples, sgpp::base::DataMatrix(dimensions, dimensions));
  for (size_t i = 0; i < samples; ++i) {
    sgpp::base::DataVector sample = distribution->operator()();
    sgpp::base::DataVector sampleGradient = gradient->gradientAt(sample);
    for (size_t d = 0; d < dimensions; ++d) {
      sgpp::base::DataVector col = sampleGradient;
      col.mult(sampleGradient[d]);
      sampleMatrices[i].setColumn(d, col);
    }
    matrix.add(sampleMatrices[i]);
  }
  matrix.mult(1.0 / static_cast<double>(samples));

  sgpp::base::DataMatrix eigenVectorMatrix(dimensions, dimensions);
  sgpp::base::DataVector eigenValues(dimensions);
  Tools::svd(matrix, eigenVectorMatrix, eigenValues);
  //size_t n = cutoff.evaluate(eigenVectorMatrix, eigenValues, sampleMatrices);
  return nullptr;
}

}  // namespace base
}  // namespace sgpp
