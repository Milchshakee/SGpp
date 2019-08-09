#pragma once

#include <random>

namespace sgpp {
namespace base {

class RandomUniformDistribution : public RandomDistribution {
 public:
  RandomUniformDistribution(size_t size, uint64_t seed, size_t dimensions);

  void generate();

 private:
  std::vector<std::uniform_real_distribution<double>> distributions;
};
}  // namespace base
}  // namespace sgpp
