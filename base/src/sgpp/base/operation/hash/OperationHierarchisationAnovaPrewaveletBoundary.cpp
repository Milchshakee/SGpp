#include <sgpp/base/operation/hash/OperationHierarchisationAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/algorithm/sweep_anova.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaLinearBoundary.hpp>

void sgpp::base::OperationHierarchisationAnovaPrewaveletBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationAnovaPrewaveletBoundary func(grid.getStorage());

  HierarchisationAnovaLinearBoundary linear(grid);
  sweep_anova<HierarchisationAnovaLinearBoundary> h(linear, grid.getStorage());
  for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
    h.sweep1D_AnovaBoundary(node_values, node_values, i);
  }

  sweep_anova<HierarchisationAnovaPrewaveletBoundary> s(func, grid.getStorage());
  for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
    s.sweep1D_AnovaBoundary(node_values, node_values, i);
  }
  DataVector copy = node_values;
}

void sgpp::base::OperationHierarchisationAnovaPrewaveletBoundary::doDehierarchisation(
    DataVector& alpha) {
  throw not_implemented_exception();
}
