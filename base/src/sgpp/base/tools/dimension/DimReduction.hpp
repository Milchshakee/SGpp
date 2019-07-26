// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DIMREDUCTION_HPP
#define DIMREDUCTION_HPP
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/optimization/function/scalar/ScalarFunction.hpp"
#include "sgpp/optimization/function/vector/VectorFunction.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {

class TransformationFunction : public sgpp::optimization::VectorFunction {
 public:
  TransformationFunction() : VectorFunction(0, 0) {};
  TransformationFunction(DataMatrix transformation);
  ~TransformationFunction() override = default;

  void eval(const base::DataVector& x, DataVector& out) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;
  size_t getOldDimensions();
  size_t getNewDimensions();

 private:
  DataMatrix transformation;
};

  class EvalFunction : public sgpp::optimization::ScalarFunction {
 public:
    EvalFunction() : ScalarFunction(0) {}
    EvalFunction(const SGridSample& sample);
  ~EvalFunction() override = default;

  double eval(const base::DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;

 private:
  const SGridSample* sample;
  };

template <class INPUT, class INFO, class OUTPUT>
class Cutter {
 public:
  virtual OUTPUT cut(const INPUT& input, const INFO& info) = 0;
};

template <class T>
class Result {
 public:
  double calcMcL2Error(optimization::ScalarFunction& func, size_t paths,
                       uint64_t seed = std::mt19937_64::default_seed)
  {
    std::mt19937_64 rand(seed);
    std::uniform_real_distribution<double> dist(0, 1);
    size_t funcDimensions = getTransformationFunction().getOldDimensions();
    size_t newDimensions = getTransformationFunction().getNewDimensions();

    sgpp::base::DataVector point(funcDimensions);
    double res = 0;

    for (size_t i = 0; i < paths; i++) {
      for (size_t d = 0; d < funcDimensions; d++) {
        point[d] = dist(rand);
      }
      double val = func.eval(point);
      DataVector out(newDimensions);
      getTransformationFunction().eval(point, out);
      res += pow(val - getReducedFunction().eval(out), 2);
    }

    return sqrt(res / static_cast<double>(paths));
  }

  virtual optimization::ScalarFunction& getReducedFunction() = 0;
  virtual TransformationFunction& getTransformationFunction() = 0;
  virtual T& getReducedOutput() = 0;
};

template <class INPUT, class INFO, class OUTPUT>
class Reducer {
 public:
  Reducer() = default;
  virtual ~Reducer() = default;

  OUTPUT evaluateAndCut(INPUT& input, Cutter<INPUT, INFO, OUTPUT>& cutter) {
    INFO i = evaluateFunction(input);
    OUTPUT result = cutter.cut(input, i);
    return result;
  }

  virtual INFO evaluate(INPUT& input) = 0;
};
}  // namespace base
}  // namespace sgpp

#endif
