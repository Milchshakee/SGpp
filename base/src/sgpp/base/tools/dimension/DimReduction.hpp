// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DIMREDUCTION_HPP
#define DIMREDUCTION_HPP
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/optimization/function/scalar/ScalarFunction.hpp"
#include "sgpp/optimization/function/vector/VectorFunction.hpp"

namespace sgpp {
namespace base {

    
class TransformationFunction : public sgpp::optimization::VectorFunction {
 public:
  TransformationFunction(DataMatrix transformation);
  ~TransformationFunction() override = default;

  void eval(const base::DataVector& x, DataVector& out) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;

 private:
  DataMatrix transformation;
};

template <class INPUT, class INFO, class OUTPUT>
class Cutter {
 public:
  virtual OUTPUT cut(const INPUT& input, const INFO& info) = 0;
};

template <class INPUT, class INFO, class OUTPUT>
class FixedCutter : public Cutter<INPUT, INFO, OUTPUT> {
 public:
  FixedCutter(size_t n) : n(n) {}


 private:
  size_t n;
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
