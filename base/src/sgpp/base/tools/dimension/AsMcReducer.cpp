#include "AsMcReducer.hpp"
#include "Tools.hpp"

namespace sgpp {
namespace base {

Sample<DataMatrix> AsMcReducer::fromGradientSample(
    const Sample<DataVector>& gradients) {
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
  return Sample<DataMatrix>(gradients.getKeys(), out);
}


Sample<DataMatrix> AsMcReducer::fromFiniteDifferences(optimization::ScalarFunction& func,
  VectorDistribution& v) {
}


AsMcIntervalCutter::AsMcIntervalCutter(size_t bootstrapSamples) : bootstrapSamples(bootstrapSamples) {
}

AsResult AsMcIntervalCutter::cut(const Sample<DataMatrix>& input, const AsInfo& info) {
  size_t dimensions = info.eigenValues.size();
  std::mt19937_64 prng;
  std::uniform_int_distribution<size_t> dist(0, input.getSize() - 1);

  std::vector<std::pair<double, double>> eigenValueIntervals(dimensions);
  for (size_t d = 0; d < dimensions; ++d) {
    eigenValueIntervals[d].first = info.eigenValues[d];
    eigenValueIntervals[d].second = info.eigenValues[d];
  }

  for (size_t i = 0; i < bootstrapSamples; i++) {
    sgpp::base::DataMatrix bootstrapMatrix(info.eigenValues.size(), info.eigenValues.size());
    for (size_t j = 0; j < input.getSize(); i++) {
      size_t l = dist(prng);
      bootstrapMatrix.add(input.getValues()[l]);
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
  return AsResult(info.eigenVectors, cutoff);
}

void AsMcReducer::evaluate(Sample<DataMatrix>& input, AsInfo& out) {
  size_t dimensions = input.getDimensions();
  sgpp::base::DataMatrix matrix(dimensions, dimensions);
  for (size_t i = 0; i < input.getSize(); ++i) {
    matrix.add(input.getValues()[i]);
  }
  matrix.mult(1.0 / static_cast<double>(input.getSize()));

  out.eigenVectors = sgpp::base::DataMatrix(dimensions, dimensions);
  out.eigenValues = sgpp::base::DataVector(dimensions);
  Tools::svd(matrix, out.eigenVectors, out.eigenValues);
}

}  // namespace base
}  // namespace sgpp
