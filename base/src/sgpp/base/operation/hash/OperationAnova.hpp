// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONVARIANCE_HPP
#define OPERATIONVARIANCE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>
#include "sgpp/datadriven/application/learnersgdeonoffparallel/AuxiliaryStructures.hpp"
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
   * Vector that holds levels for every dimension
   */
  typedef std::vector<sgpp::base::HashGridPoint::level_type> LevelVector;

  /**
   * Constructor
   */
  OperationAnova(GridStorage& gridStorage) : gridStorage(gridStorage) {}

  /**
   * Destructor
   */
  ~OperationAnova() = default;

  double calculateIncrementsVariance(const sgpp::base::DataVector& alpha,
      const std::vector<LevelVector>& levels);

  double calculateIncrementVariance(const DataVector& alpha,
                                  const LevelVector& levels);

  double calculateDimensionVariance(const DataVector& alpha,
                                    const AnovaComponent& comp);

  Sample<AnovaComponent, double> calculateAnovaOrderVariances(const DataVector& alpha);


 private:
  /// reference to the grid's GridStorage object
  GridStorage& gridStorage;

};

}  // namespace base
}  // namespace sgpp

#endif
