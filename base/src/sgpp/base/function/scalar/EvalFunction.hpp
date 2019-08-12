// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once 

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/tools/Sample.hpp>

namespace sgpp {
namespace base {

/**
 * Evaluates a hierachised grid sample at a given position.
 */
class EvalFunction : public ScalarFunction {
 public:
  EvalFunction();
  EvalFunction(const SGridSample& sample);
  EvalFunction(const SGridSample& sample, OperationEval& evalOp);
  ~EvalFunction() override = default;

  double eval(const base::DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;

 private:
  const SGridSample* sample;
  OperationEval* evalOp;
};
}  // namespace base
}  // namespace sgpp