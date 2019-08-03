// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/function/scalar/ScalarFunction.hpp"
#include "sgpp/base/function/vector/VectorFunction.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {

class MatrixFunction : public VectorFunction {
 public:
  MatrixFunction() : VectorFunction(0, 0){};
  MatrixFunction(DataMatrix transformation);
  ~MatrixFunction() override = default;

  void eval(const base::DataVector& x, DataVector& out) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;
  size_t getOldDimensions();
  size_t getNewDimensions();

 private:
  DataMatrix transformation;
};
}  // namespace base
}  // namespace sgpp