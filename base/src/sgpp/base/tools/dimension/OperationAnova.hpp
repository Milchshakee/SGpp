// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONANOVA_HPP
#define OPERATIONANOVA_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include "sgpp/base/grid/type/AnovaBoundaryGrid.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {

/**
 * Operation that cpnverts a given basis into the normal, linear hat basis and vice versa
 *
 */
class OperationAnova {
 public:

  /**
   * Constructor
   */
  OperationAnova(GridStorage& gridStorage) : gridStorage(gridStorage) {}

  /**
   * Destructor
   */
  ~OperationAnova() = default;

  Sample<AnovaBoundaryGrid::AnovaComponent, double> calculateAnovaComponentVariances(
      const DataVector& alpha);


 private:
  /// reference to the grid's GridStorage object
  GridStorage& gridStorage;

};

}  // namespace base
}  // namespace sgpp

#endif
