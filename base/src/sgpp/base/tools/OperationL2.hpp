// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>

namespace sgpp {
namespace base {

/**
 * Operation that cpnverts a given basis into the normal, linear hat basis and vice versa
 *
 */
class OperationL2 {
 public:
  /**
   * Constructor
   */
  OperationL2(uint64_t seed, size_t samples) : seed(seed), samples(samples) {}

  /**
   * Destructor
   */
  ~OperationL2() = default;

  double calculateL2Norm(ScalarFunction& func);

    double calculateMcL2Error(ScalarFunction& func, VectorFunction& transformation, ScalarFunction& reduced);

 private:
  uint64_t seed;
  size_t samples;
};

}  // namespace base
}  // namespace sgpp
