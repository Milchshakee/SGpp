// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>

namespace sgpp {
namespace base {

/**
 * Hierarchisation on ANOVA boundary grids
 *
 */
class OperationHierarchisationAnovaPrewaveletBoundary : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param grid the grid object
   */
  explicit OperationHierarchisationAnovaPrewaveletBoundary(Grid& grid) : grid(grid) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationAnovaPrewaveletBoundary() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

 protected:
   Grid& grid;

};

}  // namespace base
}  // namespace sgpp
