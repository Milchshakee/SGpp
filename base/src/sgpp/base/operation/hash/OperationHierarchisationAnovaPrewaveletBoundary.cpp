#include <sgpp/base/operation/hash/OperationHierarchisationAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/algorithm/sweep_anova.hpp>

void sgpp::base::OperationHierarchisationAnovaPrewaveletBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationAnovaPrewaveletBoundary func(grid.getStorage());
  DataVector copy = node_values;
  //func.convertLevelZeroAnsatzFunction(node_values, copy);
  sweep_anova<HierarchisationAnovaPrewaveletBoundary> s(func, grid.getStorage());
  for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
    s.sweep1D_AnovaBoundary(node_values, node_values, i);
  }
}

void sgpp::base::OperationHierarchisationAnovaPrewaveletBoundary::doDehierarchisation(
    DataVector& alpha) {
  throw not_implemented_exception();
}
