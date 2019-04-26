#include "VectorDistribution.hpp"


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

sgpp::base::DataVector UniformVectorDistribution::operator()() {
  sgpp::base::DataVector toReturn(domain.size());
  for (size_t dim = 0; dim < domain.size(); ++dim) {
    std::uniform_real_distribution<double> distribution(domain[dim].first, domain[dim].second);
    toReturn[dim] = distribution(prng);
  }
  return std::move(toReturn);
}
