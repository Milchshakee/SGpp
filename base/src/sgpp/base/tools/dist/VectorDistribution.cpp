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


RandomDistributionGenerator::RandomDistributionGenerator(size_t dimensions, std::uint64_t seed)
    : dimensions(dimensions), prng(std::mt19937_64(seed)) {}

PdfGenerator::PdfGenerator(std::uint64_t seed, ScalarFunction& pdf, size_t iterations,
  double stepSize) {
}

MultiRandomDistributionGenerator::MultiRandomDistributionGenerator(std::uint64_t seed,
  std::vector<RandomDistributionGenerator> gens) {
}

}  // namespace base
}  // namespace sgpp
