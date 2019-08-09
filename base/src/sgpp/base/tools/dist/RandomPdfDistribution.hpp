// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/dist/VectorDistribution.hpp>
#include <sgpp/base/tools/Sample.hpp>

namespace sgpp {
namespace base {

class RandomPdfDistribution : public RandomDistribution {
 public:
  RandomPdfDistribution(size_t size, size_t dimensions, uint64_t seed, ScalarFunction& pdf,
                        size_t iterations, double stepSize);

  void generate();

private:
  ScalarFunction& pdf;
  size_t iterations;
  double stepSize;

  double getProbability(DataVector& v);
  void generateVector(size_t index);
};

}  // namespace base
}  // namespace sgpp
