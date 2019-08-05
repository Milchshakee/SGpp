#include <sgpp/base/tools/dist/VectorDistribution.hpp>
#include <sgpp/base/tools/dist/RandomUniformDistribution.hpp>

namespace sgpp {
namespace base {

RandomUniformDistribution::RandomUniformDistribution(size_t size, uint64_t seed, size_t dimensions)
    : RandomDistribution(size, dimensions, seed) {
  for (size_t dim = 0; dim < dimensions; ++dim) {
    std::uniform_real_distribution<double> distribution(0, 1);
    distributions[dim] = distribution;
  }

  generate();
}

void RandomUniformDistribution::generate() {
  for (DataVector& v : vectors) {
    for (size_t dim = 0; dim < dimensions; ++dim) {
      v[dim] = distributions[dim](prng);
    }
  }
}
}  // namespace base
}  // namespace sgpp