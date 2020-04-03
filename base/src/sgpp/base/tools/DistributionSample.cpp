#include "DistributionSample.hpp"

namespace sgpp {
namespace base {


DistributionSample::DistributionSample(const std::vector<DataVector>& vectors) : size(vectors.size()), dimensions(vectors[0].getSize()), vectors(vectors) {
}

DistributionSample::DistributionSample(const DataMatrix& data)
    : size(data.getNrows()), dimensions(data.getNcols()), vectors(size) {
  for (size_t i = 0; i < size; i++) {
    data.getRow(i, vectors[i]);
  }
}

DistributionSample::DistributionSample(sgpp::base::Grid& grid) : size(grid.getSize()), dimensions(grid.getDimension()), vectors(size, DataVector(dimensions, 0.0)) {
  for (size_t i = 0; i < grid.getSize(); i++) {
    sgpp::base::GridPoint& gp = grid.getStorage().getPoint(i);
    for (size_t j = 0; j < grid.getDimension(); j++) {
      vectors[i][j] = gp.getStandardCoordinate(j);
    }
  }
}

DistributionSample::DistributionSample(size_t size, sgpp::base::DistributionsVector& dist)
    : size(size), dimensions(dist.getSize()), vectors(size) {
  for (size_t i = 0; i < size; i++) {
    vectors[i] = dist.sample();
  }
}

double getProbability(ScalarFunction& pdf, size_t dimensions, DataVector& v) {
  for (size_t d = 0; d < dimensions; ++d) {
    if (v[d] < 0.0 || v[d] > 1.0) {
      return 0.0;
    }
  }
  return pdf.eval(v);
}

DistributionSample::DistributionSample(size_t size, std::default_random_engine& gen,
                                       ScalarFunction& pdf,
  size_t iterations, double stepSize) : size(size), dimensions(pdf.getNumberOfParameters()), vectors(size) {
  std::normal_distribution<double> ndistribution(0, stepSize);
  std::uniform_real_distribution<double> udistribution(0, 1);

  for (size_t index = 0; index < size; ++index) {
    DataVector& x = vectors[index];
    for (size_t d = 0; d < dimensions; ++d) {
      x[d] = udistribution(gen);
    }

    DataVector proposed(dimensions);
    for (size_t i = 0; i < iterations; ++i) {
      double px = getProbability(pdf, dimensions, x);
      for (size_t d = 0; d < dimensions; ++d) {
        proposed[d] = x[d] + ndistribution(gen);
      }
      double pp = getProbability(pdf, dimensions, proposed);
      double alpha = 0;
      if (px == 0 && pp > 0) {
        alpha = 2.0;
      } else if (px == 0 && pp == 0) {
        alpha = 0;
      } else {
        alpha = pp / px;
      }
      if (alpha > 1 || udistribution(gen) < alpha) {
        x = proposed;
      } else {
        // Do nothing
      }
    }
  }
}

const std::vector<DataVector>& DistributionSample::getVectors() const { return vectors; }

DataMatrix DistributionSample::getAsDataMatrix() const
{
  DataMatrix data(getSize(), getDimensions());
  for (size_t i = 0; i < size; i++) {
    data.setRow(i, vectors[i]);
  }
  return data;
}

size_t DistributionSample::getSize() const { return size; }

size_t DistributionSample::getDimensions() const { return dimensions; }

}  // namespace base
}  // namespace sgpp
