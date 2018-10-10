// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>

namespace sgpp {
namespace optimization {

/**
 * Sparse grid interpolant gradient of a scalar-valued function
 * \f$f\colon \mathbb{R}^d \to \mathbb{R}\f$
 * (in contrast to ASInterpolantScalarFunctionGradient's \f$f\colon [0,1]^d \to \mathbb{R}\f$)
 *
 * @see ASInterpolantScalarFunction
 */
class ASInterpolantScalarFunctionGradient : public ScalarFunctionGradient {
 public:
  /**
   * Constructor.
   * Do not destruct the grid before the
   * ASInterpolantScalarFunctionGradient object!
   *
   * @param grid  sparse grid
   * @param alpha coefficient vector
   */
  ASInterpolantScalarFunctionGradient(base::Grid& grid, const base::DataVector& alpha)
      : ScalarFunctionGradient(grid.getDimension()),
        grid(grid),
        opEvalGradient(op_factory::createOperationEvalGradientNaive(grid)),
        alpha(alpha) {}

  /**
   * Destructor.
   */
  ~ASInterpolantScalarFunctionGradient() override {}

  /**
   * Evaluation of the function and its gradient.
   *
   * @param      x        evaluation point \f$\vec{x} \in \mathbb{R}^d\f$
   * @param[out] gradient gradient
   *                      \f$\nabla f(\vec{x}) \in \mathbb{R}^d\f$
   * @return              \f$f(\vec{x})\f$
   */
  inline double eval(const base::DataVector& x, base::DataVector& gradient) override {
    return opEvalGradient->evalGradient(alpha, x, gradient);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const override {
    clone = std::unique_ptr<ScalarFunctionGradient>(
        new ASInterpolantScalarFunctionGradient(grid, alpha));
  }

  /**
   * @return coefficient vector
   */
  const base::DataVector& getAlpha() const { return alpha; }

  /**
   * @param alpha coefficient vector
   */
  void setAlpha(const base::DataVector& alpha) { this->alpha = alpha; }

 protected:
  /// sparse grid
  base::Grid& grid;
  /// pointer to evaluation operation
  std::unique_ptr<base::OperationEvalGradient> opEvalGradient;
  /// coefficient vector
  base::DataVector alpha;
};
}  // namespace optimization
}  // namespace sgpp
