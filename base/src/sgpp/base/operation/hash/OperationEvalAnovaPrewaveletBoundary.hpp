// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>
#include <sgpp/base/grid/storage/hashmap/AnovaGridIterator.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements OperationEval for ANOVA grids with prewavelet basis ansatzfunctions with
 * boundaries
 *
 */
class OperationEvalAnovaPrewaveletBoundary : public OperationEval {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  explicit OperationEvalAnovaPrewaveletBoundary(GridStorage& storage) : storage(storage) {}

  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   * @param component the ANOVA component for which the ansatz functions should be evaluated
   */
  explicit OperationEvalAnovaPrewaveletBoundary(GridStorage& storage,
                                                AnovaBoundaryGrid::AnovaComponent& component)
      : storage(storage), component(component) {}

  /**
   * Destructor
   */
  ~OperationEvalAnovaPrewaveletBoundary() override {}

  double eval(const DataVector& alpha, const DataVector& point) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  AnovaBoundaryGrid::AnovaComponent component;
};

}  // namespace base
}  // namespace sgpp
