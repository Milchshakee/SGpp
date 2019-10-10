#include <sgpp/base/operation/hash/OperationHierarchisationAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/algorithm/sweep_anova.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationAnovaDecomposition.hpp>

void sgpp::base::OperationHierarchisationAnovaPrewaveletBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationAnovaLinearBoundary linear(grid);
  sweep_anova<HierarchisationAnovaLinearBoundary> h(linear, grid.getStorage());
  for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
    h.sweep1D_AnovaBoundary(node_values, node_values, i);
  }

  HierarchisationAnovaPrewaveletBoundary func(grid.getStorage());
  sweep_anova<HierarchisationAnovaPrewaveletBoundary> s(func, grid.getStorage());
  for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
    s.sweep1D_AnovaBoundary(node_values, node_values, i);
  }
}

void sgpp::base::OperationHierarchisationAnovaPrewaveletBoundary::doDehierarchisation(
    DataVector& alpha) {
  throw not_implemented_exception();
}
