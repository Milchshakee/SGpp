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

  std::unique_ptr<OUTPUT> evaluateAndReduce(INPUT& input, Cutter<INFO>& cutoff, INFO& out) {
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
