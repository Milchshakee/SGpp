#include "DimReduction.hpp"
#include <random>

namespace sgpp {
namespace base {

VectorDistribution::VectorDistribution(std::uint64_t seed) : prng(std::mt19937_64(seed)) {}

UniformVectorDistribution::UniformVectorDistribution(uint64_t seed, size_t dimensions)
    : UniformVectorDistribution(seed, Domain(dimensions, {0.0, 1.0})) {}

UniformVectorDistribution::UniformVectorDistribution(uint64_t seed, Domain domain)
    : VectorDistribution(seed), domain(domain), distributions(domain.size()) {
  for (size_t dim = 0; dim < domain.size(); ++dim) {
    std::uniform_real_distribution<double> distribution(domain[dim].first, domain[dim].second);
    distributions[dim] = distribution;
  }
}

DataVector UniformVectorDistribution::operator()() {
  DataVector toReturn(domain.size());
  for (size_t dim = 0; dim < domain.size(); ++dim) {
    toReturn[dim] = distributions[dim](prng);
  }
  return std::move(toReturn);
}

LengthVectorDistribution::LengthVectorDistribution(uint64_t seed, size_t dimensions, double length)
    : VectorDistribution(seed), dimensions(dimensions), length(length) {}

DataVector LengthVectorDistribution::operator()() {
  DataVector toReturn(dimensions);
  for (size_t dim = 0; dim < dimensions; ++dim) {
    toReturn[dim] = distribution(prng);
  }
  toReturn.mult(toReturn.l2Norm() / length);
  return std::move(toReturn);
}

}  // namespace base
}  // namespace sgpp
