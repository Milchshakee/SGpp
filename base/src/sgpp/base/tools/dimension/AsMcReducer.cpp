#include "AsMcReducer.hpp"
#include "EigenHelper.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {


AsMcResult::AsMcResult(const AsMcInput& input, const DataMatrix& m, size_t n)
    : AsResult<sgpp::optimization::ScalarFunction>(m, n) {
  input.function.clone(function);
}

optimization::ScalarFunction& AsMcResult::getReducedFunction() { return *function; }

optimization::ScalarFunction& AsMcResult::getReducedOutput() { return *function; }

PointSample<DataMatrix> AsMcReducer::fromGradientSample(const PointSample<DataVector>& gradients) {
  std::vector<DataMatrix> out(gradients.getSize(),
                              DataMatrix(gradients.getDimensions(), gradients.getDimensions()));
  for (size_t i = 0; i < gradients.getSize(); ++i) {
    sgpp::base::DataVector sampleGradient = gradients.getValues()[i];
    for (size_t d = 0; d < gradients.getDimensions(); ++d) {
      sgpp::base::DataVector col = sampleGradient;
      col.mult(sampleGradient[d]);
      out[i].setColumn(d, col);
    }
  }
  return PointSample<DataMatrix>(gradients.getKeys(), out);
}


PointSample<DataMatrix> AsMcReducer::fromFiniteDifferences(optimization::ScalarFunction& func,
  VectorDistribution& v) {
}

AsMcFixedCutter::AsMcFixedCutter(size_t n) : n(n) {
}

AsMcResult AsMcFixedCutter::cut(const AsMcInput& input, const AsInfo& info) {
  return AsMcResult(input, info.eigenVectors, n);
}

AsMcIntervalCutter::AsMcIntervalCutter(size_t bootstrapSamples) : bootstrapSamples(bootstrapSamples) {
}

  
AsInfo AsMcReducer::evaluate(AsMcInput& input) {
  size_t dimensions = input.samples.getDimensions();
  sgpp::base::DataMatrix matrix(dimensions, dimensions);
  for (size_t i = 0; i < input.samples.getSize(); ++i) {
    matrix.add(input.samples.getValues()[i]);
  }
  matrix.mult(1.0 / static_cast<double>(input.samples.getSize()));

  AsInfo i;
  i.eigenVectors = sgpp::base::DataMatrix(dimensions, dimensions);
  i.eigenValues = sgpp::base::DataVector(dimensions);
  i.permutation = sgpp::base::DataMatrix(dimensions, dimensions);
  EigenHelper::svd(EigenHelper::toEigen(matrix), i.eigenVectors, i.eigenValues, i.permutation);
  return i;
}

AsMcResult AsMcIntervalCutter::cut(const AsMcInput& input, const AsInfo& info) {
  size_t dimensions = info.eigenValues.size();
  std::mt19937_64 prng;
  std::uniform_int_distribution<size_t> dist(0, input.samples.getSize() - 1);

  std::vector<std::pair<double, double>> eigenValueIntervals(dimensions);
  for (size_t d = 0; d < dimensions; ++d) {
    eigenValueIntervals[d].first = info.eigenValues[d];
    eigenValueIntervals[d].second = info.eigenValues[d];
  }

  for (size_t i = 0; i < bootstrapSamples; i++) {
    sgpp::base::DataMatrix bootstrapMatrix(info.eigenValues.size(), info.eigenValues.size());
    for (size_t j = 0; j < input.samples.getSize(); i++) {
      size_t l = dist(prng);
      bootstrapMatrix.add(input.samples.getValues()[l]);
    }
    bootstrapMatrix.mult(1.0 / static_cast<double>(bootstrapSamples));

    sgpp::base::DataMatrix bootstrapEigenVectorMatrix(dimensions, dimensions);
    sgpp::base::DataVector bootstrapEigenValues(dimensions);
    sgpp::base::DataMatrix bootstrapPermutations(dimensions, dimensions);
    EigenHelper::svd(EigenHelper::toEigen(bootstrapMatrix), bootstrapEigenVectorMatrix, bootstrapEigenValues, bootstrapPermutations);
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
  return AsMcResult(input, info.eigenVectors, cutoff);
}

}  // namespace base
}  // namespace sgpp
