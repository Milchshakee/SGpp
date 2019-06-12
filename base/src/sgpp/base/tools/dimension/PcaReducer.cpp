#include <sgpp/base/tools/dimension/PcaReducer.hpp>
#include "Tools.hpp"


sgpp::base::PcaResult::PcaResult(const DataMatrix& m, size_t n) {
  transformation = m;
  transformation.resizeRowsCols(m.getNrows(), n);
}


sgpp::base::PcaFixedCutter::PcaFixedCutter(size_t n) : n(n) {
}

sgpp::base::PcaResult sgpp::base::PcaFixedCutter::cut(const VectorDistribution& input,
                                                      const PcaInfo& info){
  return PcaResult(info.eigenVectors, n);
}

sgpp::base::PcaVarianceCutter::PcaVarianceCutter(double variancePercentage) : variancePercentage(variancePercentage) {
}

sgpp::base::PcaResult sgpp::base::PcaVarianceCutter::cut(const VectorDistribution& input,
                                                         const PcaInfo& info) {
  double sum = 0;
  for (size_t d = 0; d < info.eigenValues.size(); ++d) {
    sum += info.eigenValues[d];
  }

  double percSum = 0;
  for (size_t d = 0; d < info.eigenValues.size(); ++d) {
    double varPercentage = info.eigenValues[d] / sum;
    percSum += varPercentage;
    if (percSum >= variancePercentage) {
      return PcaResult(info.eigenVectors, d + 1);
    }
  }
  return PcaResult(info.eigenVectors, info.eigenValues.size());
}

sgpp::base::FixedDistribution sgpp::base::PcaResult::apply(const VectorDistribution& input) {
  std::vector<DataVector> newDist(input.getSize());
  for (size_t c = 0; c < input.getSize(); c++) {
    newDist[c] = Tools::mult(transformation, input.getVectors()[c]);
  }
  return FixedDistribution(input.getSize(), input.getDimensions(), newDist);
}

sgpp::base::DataMatrix&& centerMean(sgpp::base::VectorDistribution& input) {
  std::vector<double> means(input.getDimensions());
  for (size_t d = 0; d < input.getDimensions(); ++d) {
    double mean = 0;
    for (size_t c = 0; c < input.getSize(); c++) {
      mean += input.getVectors()[c][d];
    }
    mean /= input.getSize();
    means[d] = mean;
  }

  sgpp::base::DataMatrix newVectors(input.getDimensions(), input.getSize());
  for (size_t c = 0; c < input.getSize(); c++) {
    newVectors.setColumn(c, input.getVectors()[c]);
    for (size_t d = 0; d < input.getDimensions(); ++d) {
      newVectors.set(d, c, newVectors.get(d, c) - means[d]);
    }
  }
  return std::move(newVectors);
}


sgpp::base::PcaReducer::PcaReducer(size_t iterations, uint64_t seed) : iterations(iterations), seed(seed) {
}

sgpp::base::PcaInfo sgpp::base::PcaReducer::evaluate(VectorDistribution& input) {
  size_t dimension = input.getDimensions();
  DataMatrix dist = centerMean(input);

  sgpp::base::DataMatrix random(dimension, dimension);
  RandomOrientationDistribution orientations(input.getDimensions(), seed, input.getDimensions(), 1);
  for (size_t d = 0; d < dimension; ++d) {
    random.setColumn(d, orientations.getVectors()[d]);
  }

  Eigen::MatrixXd x = Tools::toEigen(dist);
  Eigen::MatrixXd r = Tools::toEigen(random);
  Eigen::MatrixXd s(dimension, dimension);
  Eigen::MatrixXd e(dimension, dimension);
  for (size_t it = 0; it < iterations; ++it) {
    s = x.transpose() * (x * r);
    e = r.transpose() * s;
    r = s;
  }

  sgpp::base::DataMatrix eigenVectorMatrix = Tools::fromEigen(r);
  sgpp::base::DataVector eigenValues(dimension);
  for (size_t c = 0; c < dimension; c++) {
    eigenValues.set(c, e(c, c));
  }
  return {eigenVectorMatrix, eigenValues};
}
