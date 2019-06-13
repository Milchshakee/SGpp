// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationAnovaBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

//void OperationHierarchisationAnovaBoundary::doHierarchisation(DataVector& node_values) {
//  HierarchisationAnovaBoundary func(storage);
//  sweep<HierarchisationAnovaBoundary> s(func, storage);
//
//  // N D case
//  if (this->storage.getDimension() > 1) {
//    for (size_t i = 0; i < this->storage.getDimension(); i++) {
//      s.sweep1D_Boundary(node_values, node_values, i);
//    }
//  } else {  // 1 D case
//    s.sweep1D(node_values, node_values, 0);
//  }
//}
//
//void OperationHierarchisationAnovaBoundary::doDehierarchisation(DataVector& alpha) {
//  DehierarchisationAnovaBoundary func(storage);
//  sweep<DehierarchisationAnovaBoundary> s(func, storage);
//
//  // N D case
//  if (this->storage.getDimension() > 1) {
//    for (size_t i = 0; i < this->storage.getDimension(); i++) {
//      s.sweep1D_Boundary(alpha, alpha, i);
//    }
//  } else {  // 1 D case
//    s.sweep1D(alpha, alpha, 0);
//  }
//}

}  // namespace base
}  // namespace sgpp
