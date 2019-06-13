#include "VectorDistribution.hpp"

namespace sgpp {
namespace base {

VectorDistribution::VectorDistribution(size_t size, size_t dimensions)
    : size(size), dimensions(dimensions), vectors(size, DataVector(dimensions, 0)) {}

const std::vector<DataVector>& VectorDistribution::getVectors() const { return vectors; }

DataMatrix VectorDistribution::getAsDataMatrix() const
{
  DataMatrix data(getSize(), getDimensions());
  for (size_t i = 0; i < size; i++) {
    data.setColumn(i, vectors[i]);
  }
  return data;
}

size_t VectorDistribution::getSize() const { return size; }

size_t VectorDistribution::getDimensions() const { return dimensions; }

GridDistribution::GridDistribution(sgpp::base::Grid& grid) : VectorDistribution(grid.getSize(), grid.getDimension())
{
  for (size_t i = 0; i < grid.getSize(); i++) {
    sgpp::base::GridPoint& gp = grid.getStorage().getPoint(i);
    for (size_t j = 0; j < grid.getDimension(); j++) {
      vectors[i][j] = gp.getStandardCoordinate(j);
    }
  }
}


FixedDistribution::FixedDistribution(const DataMatrix& data)
    : VectorDistribution(data.getNcols(), data.getNrows())
{
  for (size_t i = 0; i < size; i++) {
    data.getColumn(i, vectors[i]);
  }
}

FixedDistribution::FixedDistribution(size_t size, size_t dimensions,
  const std::vector<DataVector>& vectors) : VectorDistribution(size, dimensions) {
  VectorDistribution::vectors = vectors;
}

RandomDistribution::RandomDistribution(size_t size, size_t dimensions, std::uint64_t seed)
    : VectorDistribution(size, dimensions), prng(std::mt19937_64(seed)) {}

RandomUniformDistribution::RandomUniformDistribution(size_t size, uint64_t seed, size_t dimensions)
    : RandomUniformDistribution(size, seed, Domain(dimensions, {0.0, 1.0})) {}

RandomUniformDistribution::RandomUniformDistribution(size_t size, uint64_t seed, Domain domain)
    : RandomDistribution(size, domain.size(), seed), domain(domain), distributions(domain.size()) {
  for (size_t dim = 0; dim < domain.size(); ++dim) {
    std::uniform_real_distribution<double> distribution(domain[dim].first, domain[dim].second);
    distributions[dim] = distribution;
  }

  generate();
}

void RandomUniformDistribution::generate() {
  for (DataVector& v : vectors) {
    for (size_t dim = 0; dim < domain.size(); ++dim) {
      v[dim] = distributions[dim](prng);
    }
  }
}

RandomOrientationDistribution::RandomOrientationDistribution(size_t size, uint64_t seed, size_t dimensions, double length)
    : RandomDistribution(size, dimensions, seed), dimensions(dimensions), length(length)
{
  generate();
}

void RandomOrientationDistribution::generate() {
  for (DataVector& v : vectors) {
    for (size_t dim = 0; dim < dimensions; ++dim) {
      v[dim] = distribution(prng);
    }
    v.mult(v.l2Norm() / length);
  }
}
}  // namespace base
}  // namespace sgpp
