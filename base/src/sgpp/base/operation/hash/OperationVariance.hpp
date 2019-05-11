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

namespace sgpp {
namespace base {

/**
 * Operation that cpnverts a given basis into the normal, linear hat basis and vice versa
 *
 */
class OperationVariance {
 public:
  /**
   * Vector that holds levels for every dimension
   */
  typedef std::vector<sgpp::base::HashGridPoint::level_type> LevelVector;

  /**
   * Constructor
   */
  OperationVariance() {}

  /**
   * Destructor
   */
  ~OperationVariance() {}

  void calculateIncrementVariance(sgpp::base::Grid& grid, const DataVector& alpha,
                                  const LevelVector& levels);
};

}  // namespace base
}  // namespace sgpp

#endif
