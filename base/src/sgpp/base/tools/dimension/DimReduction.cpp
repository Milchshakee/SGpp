#include "DimReduction.hpp"
#include <random>
#include "Tools.hpp"

namespace sgpp {
namespace base {

RandomDistribution::RandomDistribution(std::uint64_t seed) : prng(std::mt19937_64(seed)) {}

UniformVectorDistribution::UniformVectorDistribution(uint64_t seed, size_t dimensions)
    : UniformVectorDistribution(seed, Domain(dimensions, {0.0, 1.0})) {}

UniformVectorDistribution::UniformVectorDistribution(uint64_t seed, Domain domain)
    : RandomDistribution(seed), domain(domain), distributions(domain.size()) {
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
    : RandomDistribution(seed), dimensions(dimensions), length(length) {}

DataVector LengthVectorDistribution::operator()() {
  DataVector toReturn(dimensions);
  for (size_t dim = 0; dim < dimensions; ++dim) {
    toReturn[dim] = distribution(prng);
  }
  toReturn.mult(toReturn.l2Norm() / length);
  return std::move(toReturn);
}

ReducedFunction::ReducedFunction(std::unique_ptr<sgpp::optimization::ScalarFunction>&& function,
  DataMatrix transformation) : ScalarFunction(transformation.getNcols()), function(std::move(function)), transformation(transformation) {}

double ReducedFunction::eval(const base::DataVector& x) {
  DataVector y = Tools::mult(transformation, x);
  return function->eval(y);
}

void ReducedFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  std::unique_ptr<ScalarFunction> ptr;
  function->clone(ptr);
  clone = std::make_unique<ReducedFunction>(std::move(ptr), transformation);
}

}  // namespace base
}  // namespace sgpp
