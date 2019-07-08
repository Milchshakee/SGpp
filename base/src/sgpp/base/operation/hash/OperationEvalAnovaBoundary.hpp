// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALANOVABOUNDARY_HPP
#define OPERATIONEVALANOVABOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements OperationEval for a grids with linear basis ansatzfunctions with
 * boundaries
 *
 */
class OperationEvalAnovaBoundary : public OperationEval {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   */
  explicit OperationEvalAnovaBoundary(GridStorage& storage) : storage(storage) {}

  /**
   * Destructor
   */
  ~OperationEvalAnovaBoundary() override {}

  double eval(const DataVector& alpha, const DataVector& point) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif
