// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HIERARCHISATIONANOVABOUNDARY_HPP
#define HIERARCHISATIONANOVABOUNDARY_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/dimension/AnovaHelper.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {
/**
 * Class that implements the hierarchisation on a linear sparse grid with boundaries. Therefore
 * the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal
 */
class HierarchisationAnovaBoundary {
 protected:
  typedef AnovaHelper::AnovaGridIterator grid_iterator;

 public:
  /**
   * Constructor, must be bind to a grid
   *
   * @param storage the grid storage object of the the grid, on which the hierarchisation should be
   * executed
   */
  explicit HierarchisationAnovaBoundary(GridStorage& storage);

  /**
   * Destructor
   */
  ~HierarchisationAnovaBoundary();

  /**
   * Implements operator() needed by the sweep class during the grid traversal. This function
   * is applied to the whole grid.
   *
   * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
   * If one is missing this code might produce some bad errors (segmentation fault, wrong
   * calculation result) So please assure that both functions do exist!
   *
   * @param source this DataVector holds the node base coefficients of the function that should be
   * applied to the sparse grid
   * @param result this DataVector holds the linear base coefficients of the sparse grid's
   * ansatz-functions
   * @param index a iterator object of the grid
   * @param dim current fixed dimension of the 'execution direction'
   */
  virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index,
                          size_t dim);

  virtual void operator()(DataVector& source, DataVector& result, HashGridIterator& index, size_t dim);

    /**
   * Recursive hierarchisaton algorithm, this algorithms works in-place -> source should be equal to
   * result
   *
   * @param source this DataVector holds the node base coefficients of the function that should be
   * applied to the sparse grid
   * @param result this DataVector holds the linear base coefficients of the sparse grid's
   * ansatz-functions
   * @param index a iterator object of the grid
   * @param dim current fixed dimension of the 'execution direction'
   * @param fl left value of the current region regarded in this step of the recursion
   * @param fr right value of the current region regarded in this step of the recursion
   */
  void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl,
           double fr, double constant);

  void hierarchiseConstant(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);

  void hierarchiseConstantRec(DataVector& source, DataVector& result, grid_iterator& index,
                              size_t dim, double constant);

private:
  GridStorage& storage;
};

}  // namespace base
}  // namespace sgpp

#endif
