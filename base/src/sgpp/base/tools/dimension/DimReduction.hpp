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

template <class T>
class CutoffCriterion {
 public:
  virtual size_t evaluate(const T& info) = 0;
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

template <class T>
class Reducer {
 public:
  Reducer(std::shared_ptr<CutoffCriterion<T>>& cutoff) : cutoff(cutoff) {}
  virtual ~Reducer() = default;

 protected:
  std::shared_ptr<CutoffCriterion<T>> cutoff;
};

class DataReducer {
 public:
  virtual std::unique_ptr<sgpp::optimization::ScalarFunction> reduceData(
      sgpp::base::DataMatrix& input) = 0;
  std::unique_ptr<sgpp::optimization::ScalarFunction> reduceFunction(
      sgpp::optimization::ScalarFunction& input, VectorDistribution& dist, size_t samples);
};

template <class T>
class FunctionReducer : public Reducer<T> {
 public:
  FunctionReducer(std::shared_ptr<CutoffCriterion<T>>& cutoff) : Reducer<T>(cutoff) {}

  std::unique_ptr<sgpp::optimization::ScalarFunction> reduceFunction(
      sgpp::optimization::ScalarFunction& input, T& out) {
    evaluateFunction(input, out);
    size_t n = Reducer<T>::cutoff->evaluate(out);
    return reduce(input, n, out);
  }

 protected:
  virtual void evaluateFunction(sgpp::optimization::ScalarFunction& input, T& out) = 0;
  virtual std::unique_ptr<sgpp::optimization::ScalarFunction> reduce(
      sgpp::optimization::ScalarFunction& input, size_t n, const T& info) = 0;
};
}  // namespace base
}  // namespace sgpp

#endif
