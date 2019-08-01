// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationQuadratureAnovaBoundary.hpp>

#include <sgpp/globaldef.hpp>
#include "sgpp/base/grid/type/AnovaBoundaryGrid.hpp"

namespace sgpp {
namespace base {

double OperationQuadratureAnovaBoundary::doQuadrature(DataVector& alpha) {
  double res = 0;
  double tmp = 0;
  GridStorage::point_type index;
  GridStorage::grid_map_iterator end_iter = storage.end();

  for (GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) {
    size_t levelSum = 0;
    size_t level0Count = 0;
    for (size_t d = 0; d < iter->first->getDimension(); d++) {
      levelSum += std::max<size_t>(iter->first->getLevel(d), 1) - 1;
      if (sgpp::base::AnovaBoundaryGrid::fromNormalLevel(iter->first->getLevel(d)) == 0) {
        level0Count++;
      }
    }

    tmp = pow(2.0, -static_cast<double>(levelSum)) * alpha.get(iter->second);
    tmp *= (pow(2.0, -static_cast<double>(level0Count)));
    res += tmp;
  }

  // multiply with determinant of "unit cube -> BoundingBox" transformation
  for (size_t d = 0; d < storage.getDimension(); d++) {
    res *= storage.getBoundingBox()->getIntervalWidth(d);
  }

  return res;
}

}  // namespace base
}  // namespace sgpp
