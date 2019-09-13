// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationAnovaLinearBoundary.hpp>
#include <sgpp/base/algorithm/sweep_anova.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

namespace sgpp {
namespace base {

void OperationHierarchisationAnovaLinearBoundary::doHierarchisation(DataVector& node_values) {}

void OperationHierarchisationAnovaLinearBoundary::doDehierarchisation(DataVector& alpha) {
  throw not_implemented_exception();
}

}  // namespace base
}  // namespace sgpp
