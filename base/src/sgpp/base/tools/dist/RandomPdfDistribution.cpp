#include <sgpp/base/tools/dist/RandomPdfDistribution.hpp>

sgpp::base::RandomPdfDistribution::RandomPdfDistribution(size_t size, size_t dimensions,
                                                         uint64_t seed, ScalarFunction& pdf,
                                                         size_t iterations, double stepSize)
    : RandomDistribution(size, dimensions, seed),
      pdf(pdf),
      iterations(iterations),
      stepSize(stepSize) {
  generate();
}

void sgpp::base::RandomPdfDistribution::generate() {
  for (size_t i = 0; i < size; ++i) {
    generateVector(i);
  }
}

void sgpp::base::RandomPdfDistribution::generateVector(size_t index) {
  std::normal_distribution<double> ndistribution(0, stepSize);
  std::uniform_real_distribution<double> udistribution(0, 1);

  DataVector& x = vectors[index];
  for (size_t d = 0; d < dimensions; ++d) {
    x[d] = udistribution(prng);
  }

  DataVector proposed(dimensions);
  for (size_t i = 0; i < iterations; ++i) {
    double px = getProbability(x);
    for (size_t d = 0; d < dimensions; ++d) {
      proposed[d] = x[d] + ndistribution(prng);
    }
    double pp = getProbability(proposed);
    double alpha = 0;
    if (px == 0 && pp > 0) {
      alpha = 2.0;
    } else if (px == 0 && pp == 0) {
      alpha = 0;
    } else {
      alpha = pp / px;
    }
    if (alpha > 1 || udistribution(prng) < alpha) {
      x = proposed;
    } else {
      // Do nothing
    }
  }
}

double sgpp::base::RandomPdfDistribution::getProbability(DataVector& v) {
  for (size_t d = 0; d < dimensions; ++d) {
    if (v[d] < 0.0 || v[d] > 1.0) {
      return 0.0;
    }
  }
  return pdf.eval(v);
}
