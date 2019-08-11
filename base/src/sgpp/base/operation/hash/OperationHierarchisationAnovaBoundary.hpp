// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONANOVABOUNDARY_HPP
#define OPERATIONHIERARCHISATIONANOVABOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>

namespace sgpp {
namespace base {

/**
 * Hierarchisation on ANOVA boundary grids
 *
 */
class OperationHierarchisationAnovaBoundary : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param grid the grid object
   */
  explicit OperationHierarchisationAnovaBoundary(Grid& grid) : grid(grid) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationAnovaBoundary() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  Grid& grid;

 public:
  /**
   * Creates an anchor for the point at level -1 in all dimensions.
   */
  static std::vector<AnovaBoundaryGrid::LevelIndexPair> defaultAnchor(size_t dims) {
    return std::vector<AnovaBoundaryGrid::LevelIndexPair>(dims, {-1, 0});
  }

  /**
   * Make all hierarchisation operations use the point at level -1 in all dimensions as an anchor.
   */
  static void setDefaultPolicy() {
    integral = false;
    anchor.clear();
  }

  /**
   * Make all hierarchisation operations set the value of the constant function to the expected value of the function (the integral of the function).
   */
  static void setIntegralPolicy() {
    integral = true;
    anchor.clear();
  }

    /**
   * Make all hierarchisation operations use a custom anchor point.
   */
  static void setAnchorPolicy(std::vector<AnovaBoundaryGrid::LevelIndexPair>& anchor) {
    OperationHierarchisationAnovaBoundary::anchor = anchor;
  }

 private:
  double getAnchorValue(std::vector<AnovaBoundaryGrid::LevelIndexPair>& anchor,
                        DataVector& node_values);

  static bool integral;
  static std::vector<AnovaBoundaryGrid::LevelIndexPair> anchor;
};

}  // namespace base
}  // namespace sgpp

#endif
