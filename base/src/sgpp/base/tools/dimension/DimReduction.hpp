// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DIMREDUCTION_HPP
#define DIMREDUCTION_HPP

#include "sgpp/optimization/function/scalar/ScalarFunction.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include <vector>
#include <random>

namespace sgpp {
namespace base {

typedef std::vector<std::pair<double, double>> Domain;

class VectorDistribution {
 public:
  VectorDistribution(std::uint64_t seed);

  virtual DataVector operator()() = 0;

 protected:
  std::mt19937_64 prng;
};

class UniformVectorDistribution : public VectorDistribution {
 public:
  UniformVectorDistribution(uint64_t seed, size_t dimensions);
  UniformVectorDistribution(uint64_t seed, Domain domain);

  DataVector operator()() override;

 private:
  Domain domain;
  std::vector<std::uniform_real_distribution<double>> distributions;
};

class LengthVectorDistribution : public VectorDistribution {
 public:
  LengthVectorDistribution(uint64_t seed, size_t dimensions, double length);

  DataVector operator()() override;

 private:
  size_t dimensions;
  double length;
  std::uniform_real_distribution<double> distribution;
};

  class ReducedFunction : public sgpp::optimization::ScalarFunction
  {
 public:
    ReducedFunction(std::unique_ptr<sgpp::optimization::ScalarFunction>&& function, DataMatrix transformation);
    ~ReducedFunction() override;

    double eval(const base::DataVector& x) override;
    void clone(std::unique_ptr<ScalarFunction>& clone) const override;
  private:
    std::unique_ptr<sgpp::optimization::ScalarFunction> function;
   DataMatrix transformation;
  };

class DataReducer {
 public:
  virtual std::unique_ptr<sgpp::optimization::ScalarFunction> reduceData(
      sgpp::base::DataMatrix& input) = 0;
  std::unique_ptr<sgpp::optimization::ScalarFunction> reduceFunction(
      sgpp::optimization::ScalarFunction& input, VectorDistribution& dist, size_t samples);
};

class FunctionReducer
{
public:
  virtual std::unique_ptr<sgpp::optimization::ScalarFunction> reduceFunction(sgpp::optimization::ScalarFunction& input) = 0;
};

}  // namespace base
}  // namespace sgpp

#endif
