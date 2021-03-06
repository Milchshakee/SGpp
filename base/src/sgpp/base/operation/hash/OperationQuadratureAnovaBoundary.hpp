// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on ANOVA boundary grids
 */
class OperationQuadratureAnovaBoundary : public OperationQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureLinear
   *
   * @param storage Pointer to the grid's GridStorage object
   */
  explicit OperationQuadratureAnovaBoundary(GridStorage& storage) : storage(storage) {}

  ~OperationQuadratureAnovaBoundary() override {}

  /**
   * Quadrature for piecewise linear hat basis functions.
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp
