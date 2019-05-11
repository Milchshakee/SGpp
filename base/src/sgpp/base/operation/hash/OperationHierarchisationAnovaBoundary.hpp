// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONHIERARCHISATIONANOVABOUNDARY_HPP
#define OPERATIONHIERARCHISATIONANOVABOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Hierarchisation on sparse grid, Anova case with boundaries
 *
 */
class OperationHierarchisationAnovaBoundary : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  explicit OperationHierarchisationAnovaBoundary(GridStorage& storage) : storage(storage) {}

  /**
   * Destructor
   */
  ~OperationHierarchisationAnovaBoundary() override {}

  void doHierarchisation(DataVector& node_values) override;
  void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif
