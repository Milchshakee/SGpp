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
    data.setRow(i, vectors[i]);
  }
  return data;
}

size_t VectorDistribution::getSize() const { return size; }

size_t VectorDistribution::getDimensions() const { return dimensions; }

RandomDistribution::RandomDistribution(size_t size, size_t dimensions, std::uint64_t seed)
    : VectorDistribution(size, dimensions), prng(std::mt19937_64(seed)) {}
}  // namespace base
}  // namespace sgpp
