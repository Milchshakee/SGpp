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
  return Sample<DataMatrix>(gradients.getVectors(), out);
}


Sample<DataMatrix> AsMcReducer::fromFiniteDifferences(optimization::ScalarFunction& func,
  VectorDistribution& v) {
}

AsMcReducer::IntervalCutoff::IntervalCutoff(size_t bootstrapSamples) : bootstrapSamples(bootstrapSamples) {}

size_t AsMcReducer::IntervalCutoff::evaluate(const McActiveSubspaceInfo& info) {
  size_t dimensions = info.eigenValues.size();
  std::mt19937_64 prng;
  std::uniform_int_distribution<size_t> dist(0, info.sampleMatrices.getSize() - 1);

  std::vector<std::pair<double, double>> eigenValueIntervals(dimensions);
  for (size_t d = 0; d < dimensions; ++d) {
    eigenValueIntervals[d].first = info.eigenValues[d];
    eigenValueIntervals[d].second = info.eigenValues[d];
  }

  for (size_t i = 0; i < bootstrapSamples; i++) {
    sgpp::base::DataMatrix bootstrapMatrix(info.eigenValues.size(), info.eigenValues.size());
    for (size_t j = 0; j < info.sampleMatrices.getSize(); i++) {
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

void AsMcReducer::evaluate(
    Sample<DataMatrix>& input, McActiveSubspaceInfo& out) {
  size_t dimensions = input.getDimensions();
  sgpp::base::DataMatrix matrix(dimensions, dimensions);
  out.sampleMatrices = input;
  for (size_t i = 0; i < input.getSize(); ++i) {
    matrix.add(out.sampleMatrices.getValues()[i]);
  }
  matrix.mult(1.0 / static_cast<double>(input.getSize()));

  out.eigenVectors = sgpp::base::DataMatrix(dimensions, dimensions);
  out.eigenValues = sgpp::base::DataVector(dimensions);
  Tools::svd(matrix, out.eigenVectors, out.eigenValues);
}

    ActiveSubspaceResult AsMcReducer::reduce(
        Sample<DataMatrix>& input, size_t c, const McActiveSubspaceInfo& info) {
  DataMatrix m = info.eigenVectors;
  m.resizeRowsCols(input.getDimensions(), c);
  return ActiveSubspaceResult {m};
    }

}  // namespace base
}  // namespace sgpp
