// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationAnovaBoundary.hpp>

#include <sgpp/globaldef.hpp>
#include "common/algorithm_sweep/HierarchisationAnovaBoundary.hpp"
#include "sgpp/base/algorithm/sweep_anova.hpp"

namespace sgpp {
namespace base {

bool OperationHierarchisationAnovaBoundary::integral = false;
std::vector<AnovaBoundaryGrid::LevelIndexPair> OperationHierarchisationAnovaBoundary::anchor;

double OperationHierarchisationAnovaBoundary::getAnchorValue(
    std::vector<AnovaBoundaryGrid::LevelIndexPair>& anchor, DataVector& node_values) {
  for (size_t i = 0; i < grid.getStorage().getSize(); ++i) {
    GridPoint& gp = grid.getStorage().getPoint(i);
    for (size_t d = 0; d < grid.getStorage().getDimension(); ++d) {
      AnovaBoundaryGrid::level_t l;
      index_t i;
      AnovaBoundaryGrid::fromNormalGridPointLevelIndex(gp.getLevel(d), gp.getIndex(d), l, i);
      if (l != anchor[d].level || i != anchor[d].index) {
        break;
      }
      return node_values[i];
    }
  }
  throw std::invalid_argument("Can't find anchor value");
}

void OperationHierarchisationAnovaBoundary::doHierarchisation(DataVector& node_values) {
  HierarchisationAnovaBoundary func(grid);
  sweep_anova<HierarchisationAnovaBoundary> s(func, grid.getStorage());
  for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
    s.sweep1D_AnovaBoundary(node_values, node_values, i);
  }
  double offset = 0;
  if (integral) {
    auto a = defaultAnchor(grid.getDimension());
    offset = getAnchorValue(a, node_values);

  } else if (!anchor.empty()) {
    offset = getAnchorValue(anchor, node_values);
  } else {
    offset = node_values[0];
  }
  double change = offset - node_values[0];
  for (size_t i = 0; i < node_values.size(); ++i) {
    node_values[i] += change;
  }
}

void OperationHierarchisationAnovaBoundary::doDehierarchisation(DataVector& alpha) {
  // DehierarchisationAnovaBoundary func(storage);
  // sweep<DehierarchisationAnovaBoundary> s(func, storage);

  //// N D case
  // if (this->storage.getDimension() > 1) {
  //  for (size_t i = 0; i < this->storage.getDimension(); i++) {
  //    s.sweep1D_Boundary(alpha, alpha, i);
  //  }
  //} else {  // 1 D case
  //  s.sweep1D(alpha, alpha, 0);
  //}
}

}  // namespace base
}  // namespace sgpp
