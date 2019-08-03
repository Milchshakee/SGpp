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

class EvalFunction : public ScalarFunction {
 public:
  EvalFunction() : ScalarFunction(0) {}
  EvalFunction(const SGridSample& sample);
  ~EvalFunction() override = default;

  double eval(const base::DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;

 private:
  const SGridSample* sample;
};
}  // namespace base
}  // namespace sgpp