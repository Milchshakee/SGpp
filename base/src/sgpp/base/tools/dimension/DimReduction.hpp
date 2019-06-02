// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DIMREDUCTION_HPP
#define DIMREDUCTION_HPP

#include <random>
#include <vector>
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/optimization/function/scalar/ScalarFunction.hpp"

namespace sgpp {
namespace base {

typedef std::vector<std::pair<double, double>> Domain;

class RandomDistribution {
 public:
  RandomDistribution(std::uint64_t seed);

  virtual DataVector operator()() = 0;

 protected:
  std::mt19937_64 prng;
};

class UniformVectorDistribution : public RandomDistribution {
 public:
  UniformVectorDistribution(uint64_t seed, size_t dimensions);
  UniformVectorDistribution(uint64_t seed, Domain domain);

  DataVector operator()() override;

 private:
  Domain domain;
  std::vector<std::uniform_real_distribution<double>> distributions;
};

class LengthVectorDistribution : public RandomDistribution {
 public:
  LengthVectorDistribution(uint64_t seed, size_t dimensions, double length);

  DataVector operator()() override;

 private:
  size_t dimensions;
  double length;
  std::uniform_real_distribution<double> distribution;
};

template <class I>
class CutoffCriterion {
 public:
  virtual size_t evaluate(const I& info) = 0;
};

template <class T>
class FixedCutoff : public CutoffCriterion<T> {
 public:
  FixedCutoff(size_t n) : n(n) {}

            size_t evaluate(const T& info) override {
    return n;
  }

 private:
  size_t n;
};

class ReducedFunction : public sgpp::optimization::ScalarFunction {
 public:
  ReducedFunction(std::unique_ptr<sgpp::optimization::ScalarFunction>&& function,
                  DataMatrix transformation);
  ~ReducedFunction() override = default;

  double eval(const base::DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;

 private:
  std::unique_ptr<sgpp::optimization::ScalarFunction> function;
  DataMatrix transformation;
};

template <class INPUT, class INFO, class OUTPUT>
class Reducer {
 public:
  Reducer() = default;
  virtual ~Reducer() = default;

  std::unique_ptr<OUTPUT> evaluateAndReduce(INPUT& input, CutoffCriterion<INFO>& cutoff, INFO& out) {
    evaluateFunction(input, out);
    size_t c = cutoff.evaluate(out);
    return reduce(input, c, out);
  }

  virtual void evaluate(INPUT& input, INFO& out) = 0;
  virtual OUTPUT reduce(
      INPUT& input, size_t c, const INFO& info) = 0;

};
}  // namespace base
}  // namespace sgpp

#endif
