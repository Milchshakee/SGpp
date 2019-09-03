// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationAnovaBoundary.hpp>
#include <sgpp/base/algorithm/sweep_anova.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaBoundary.hpp>

namespace sgpp {
namespace base {

bool OperationHierarchisationAnovaBoundary::integral = false;
std::vector<AnovaTypes::LevelIndexPair> OperationHierarchisationAnovaBoundary::anchor;

double OperationHierarchisationAnovaBoundary::getAnchorValue(
    std::vector<AnovaTypes::LevelIndexPair>& anchor, DataVector& node_values) {
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
  if (!anchor.empty()) {
    double a = getAnchorValue(anchor, node_values);
    node_values[0] = a;
  }

  DataVector copy;
  if (integral) {
    copy = node_values;
    }

  HierarchisationAnovaBoundary func(grid);
  sweep_anova<HierarchisationAnovaBoundary> s(func, grid.getStorage());
  for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
    s.sweep1D_AnovaBoundary(node_values, node_values, i);
  }

  if (integral) {
    std::unique_ptr<OperationQuadrature> opQ(
        op_factory::createOperationQuadrature(grid));
    double q = opQ->doQuadrature(node_values);
    node_values[0] = q;

    HierarchisationAnovaBoundary func2(grid);
    sweep_anova<HierarchisationAnovaBoundary> s2(func2, grid.getStorage());
    for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
      s2.sweep1D_AnovaBoundary(copy, copy, i);
    }
    node_values = copy;
  }

  if (grid.getType() == GridType::AnovaPrewaveletBoundary)
  {
    ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary func3(grid.getStorage());
    func3.convertLevelZeroAnsatzFunction(node_values, node_values);
    sweep_anova<ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary> s3(func3, grid.getStorage());
    for (size_t i = 0; i < grid.getStorage().getDimension(); i++) {
      s3.sweep1D_AnovaBoundary(node_values, node_values, i);
    }
  } }

void OperationHierarchisationAnovaBoundary::doDehierarchisation(DataVector& alpha) {
  throw not_implemented_exception();
}

}  // namespace base
}  // namespace sgpp
