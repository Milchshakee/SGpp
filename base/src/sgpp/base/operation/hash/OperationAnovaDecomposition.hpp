// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>

namespace sgpp {
namespace base {

/**
 * FirstMomemnt of sparse grid function, linear grid without boundaries
 */
class OperationAnovaDecomposition {
 public:
  /**
   * Constructor of OperationFirstMomentLinearBoundary
   *
   * @param storage Pointer to the grid's GridStorage object
   */
  explicit OperationAnovaDecomposition(GridStorage& storage) : storage(storage) {}

  ~OperationAnovaDecomposition() {}

  
  void decompose(const DataVector& alpha, DataVector& result);

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp
