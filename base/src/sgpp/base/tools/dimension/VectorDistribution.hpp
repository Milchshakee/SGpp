// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef VECTORDISTRIBUTION_HPP
#define VECTORDISTRIBUTION_HPP

#include <random>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>

typedef std::vector<std::pair<double, double>> Domain;

class VectorDistribution {
public:
  VectorDistribution(std::uint64_t seed);
  virtual ~VectorDistribution();

  virtual sgpp::base::DataVector operator()();

protected:
  std::mt19937_64 prng;
};

class UniformVectorDistribution : VectorDistribution {
 public:
  UniformVectorDistribution(uint64_t seed, size_t dimensions);
  UniformVectorDistribution(uint64_t seed, Domain domain);

  sgpp::base::DataVector operator()() override;

 private:
  Domain domain;
  std::vector<std::uniform_real_distribution<double>> distributions;
};

#endif