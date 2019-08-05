// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>

namespace sgpp {
namespace base {

  /**
   * Multiplies a given matrix with the input vector and returns the result.
   */
class MatrixFunction : public VectorFunction {
 public:
  MatrixFunction();
  MatrixFunction(DataMatrix transformation);
  ~MatrixFunction() override = default;

  void eval(const base::DataVector& x, DataVector& out) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;

 private:
  DataMatrix transformation;
};
}  // namespace base
}  // namespace sgpp