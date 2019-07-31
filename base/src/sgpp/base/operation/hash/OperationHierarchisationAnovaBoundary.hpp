// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONANOVABOUNDARY_HPP
#define OPERATIONHIERARCHISATIONANOVABOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>

#include <sgpp/globaldef.hpp>
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/grid/type/AnovaBoundaryGrid.hpp"

namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, Anova case with boundaries
 *
 */
class OperationHierarchisationAnovaBoundary : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
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
  static std::vector<AnovaBoundaryGrid::LevelIndexPair> defaultAnchor(size_t dims) {
    return std::vector<AnovaBoundaryGrid::LevelIndexPair>(dims, {-1, 0});
  }

  static void setDefaultPolicy() {
    integral = false;
    anchor.clear();
  }

  static void setIntegralPolicy() {
    integral = true;
    anchor.clear();
  }

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
