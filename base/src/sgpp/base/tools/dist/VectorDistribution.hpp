// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <random>
#include <vector>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>

namespace sgpp {
namespace base {

class VectorDistribution {
 public:
  VectorDistribution(size_t size, size_t dimensions);

  size_t getSize() const;
  size_t getDimensions() const;
  const std::vector<DataVector>& getVectors() const;
  DataMatrix getAsDataMatrix() const;

 protected:
  size_t size;
  size_t dimensions;
  std::vector<DataVector> vectors;
};

class RandomDistribution : public VectorDistribution {
 public:
  RandomDistribution(size_t size, size_t dimensions, std::uint64_t seed);

  virtual void generate() = 0;

 protected:
  std::mt19937_64 prng;
};
}  // namespace base
}  // namespace sgpp

