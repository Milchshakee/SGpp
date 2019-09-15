// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationAnovaLinearBoundary.hpp>
#include <sgpp/base/algorithm/sweep_anova.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaLinearBoundary.hpp>

namespace sgpp {
namespace base {

void OperationHierarchisationAnovaLinearBoundary::doHierarchisation(DataVector& node_values) {
  HierarchisationAnovaLinearBoundary linear(grid);
  sweep_anova<HierarchisationAnovaLinearBoundary> h(linear, grid.getStorage());
  for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
    h.sweep1D_AnovaBoundary(node_values, node_values, i);
  }
}

void OperationHierarchisationAnovaLinearBoundary::doDehierarchisation(DataVector& alpha) {
  throw not_implemented_exception();
}

}  // namespace base
}  // namespace sgpp
