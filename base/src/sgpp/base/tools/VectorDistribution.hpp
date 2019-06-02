// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef VECTORDISTRIBUTION_HPP
#define VECTORDISTRIBUTION_HPP

#include <random>
#include <vector>
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/optimization/function/scalar/ScalarFunction.hpp"
#include "sgpp/base/grid/Grid.hpp"

namespace sgpp {
namespace base {

class VectorDistribution {
 public:
  VectorDistribution(size_t size, size_t dimensions);

  size_t getSize();
  size_t getDimensions();
  const std::vector<DataVector>& getVectors();

 protected:
  size_t size;
  size_t dimensions;
  std::vector<DataVector> vectors;
};

    class GridDistribution : public VectorDistribution {
 public:
      GridDistribution(sgpp::base::Grid& grid);
    };

  class FixedDistribution : public VectorDistribution {
 public:
    FixedDistribution(size_t size, size_t dimensions, const std::vector<DataVector>& vectors);
};

class RandomDistribution : public VectorDistribution {
 public:
  RandomDistribution(size_t size, size_t dimensions, std::uint64_t seed);

  virtual void generate() = 0;

 protected:
  std::mt19937_64 prng;
};

class RandomUniformDistribution : public RandomDistribution {
 public:
  typedef std::vector<std::pair<double, double>> Domain;

  RandomUniformDistribution(size_t size, uint64_t seed, size_t dimensions);
  RandomUniformDistribution(size_t size, uint64_t seed, Domain domain);

  void generate();

 private:
  Domain domain;
  std::vector<std::uniform_real_distribution<double>> distributions;
};

class RandomOrientationDistribution : public RandomDistribution {
 public:
  RandomOrientationDistribution(size_t size, uint64_t seed, size_t dimensions, double length);

  void generate();

 private:
  size_t dimensions;
  double length;
  std::uniform_real_distribution<double> distribution;
};

}  // namespace base
}  // namespace sgpp

#endif
