// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HIERARCHISATIONANOVABOUNDARY_HPP
#define HIERARCHISATIONANOVABOUNDARY_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>
#include <sgpp/base/grid/storage/hashmap/AnovaGridIterator.hpp>

namespace sgpp {
namespace base {
/**
 * Class that implements the hierarchisation on a ANOVA grid with boundaries. Therefore
 * the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal.
 */
class HierarchisationAnovaBoundary {
 protected:
  /// Custom iterator
  typedef AnovaGridIterator grid_iterator;

 public:
  /**
   * Constructor, must be bind to a grid
   *
   * @param grid the grid object, on which the hierarchisation should be
   * executed
   */
  explicit HierarchisationAnovaBoundary(Grid& grid);


  /**
   * Destructor
   */
  ~HierarchisationAnovaBoundary();

  /**
   * Implements operator() needed by the sweep class during the grid traversal. This function
   * is applied to the whole grid.
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
           double fr);

  void hierarchiseConstant(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);

  void hierarchiseConstantRec(DataVector& source, DataVector& result, grid_iterator& index,
                              size_t dim, double constant);

private:
  Grid& grid;
};

}  // namespace base
}  // namespace sgpp

#endif
