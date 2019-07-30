// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SWEEPANOVA_HPP
#define SWEEPANOVA_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/dimension/AnovaHelper.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <utility>
#include <vector>

namespace sgpp {
namespace base {

/**
 * Standard sweep operation
 * FUNC should be a class with overwritten operator(). For an example see laplace_up_functor in
 * laplace.hpp. It must be default constructable or copyable. STORAGE must provide a grid_iterator
 * supporting left_child, step_right, up, hint and seq.
 */
template <class FUNC>
class sweep_anova {
 protected:
  typedef AnovaHelper::AnovaGridIterator grid_iterator;

  /// Object of FUNC, this is executed by sweep
  FUNC functor;
  /// Pointer to the grid's storage object
  GridStorage& storage;
  /// algorithmic dimensions, operator is applied in this dimensions
  const std::vector<size_t> algoDims;
  /// number of algorithmic dimensions
  const size_t numAlgoDims_;

 public:
  /**
   * Create a new sweep object with a default constructed functor
   *
   * @param storage the storage that contains the grid points
   */
  explicit sweep_anova(GridStorage& storage)
      : functor(),
        storage(storage),
        algoDims(storage.getAlgorithmicDimensions()),
        numAlgoDims_(storage.getAlgorithmicDimensions().size()) {}

  /**
   * Create a new sweep object with a copied functor
   *
   * @param functor the functor that is executed on the grid
   * @param storage the storage that contains the grid points
   */
  sweep_anova(FUNC& functor, GridStorage& storage)
      : functor(functor),
        storage(storage),
        algoDims(storage.getAlgorithmicDimensions()),
        numAlgoDims_(storage.getAlgorithmicDimensions().size()) {}

  /**
   * Destructor
   */
  ~sweep_anova() {}

  /**
   * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep
   * Boundaries are regarded
   *
   * @param source a DataVector containing the source coefficients of the grid points
   * @param result a DataVector containing the result coefficients of the grid points
   * @param dim_sweep the dimension in which the functor is executed
   */
  void sweep1D_AnovaBoundary(DataVector& source, DataVector& result, size_t dim_sweep) {
    // generate a list of all dimension (-dim_sweep) from
    // dimension recursion unrolling
    std::vector<size_t> dim_list;

    for (size_t i = 0; i < storage.getDimension(); i++) {
      if (i != dim_sweep) {
        dim_list.push_back(i);
      }
    }

    grid_iterator index(storage);
    index.resetToLevelMinusOne();

    sweep_AnovaBoundary_rec(source, result, index, dim_list, storage.getDimension() - 1, dim_sweep);
  }

 protected:
  /**
   * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep.
   * Boundaries are regarded
   *
   * @param source coefficients of the sparse grid
   * @param result coefficients of the function computed by sweep
   * @param index current grid position
   * @param dim_list list of dimensions, that should be handled
   * @param dim_rem number of remaining dims
   * @param dim_sweep static dimension, in this dimension the functor is executed
   */
  void sweep_AnovaBoundary_rec(DataVector& source, DataVector& result,
                               AnovaHelper::AnovaGridIterator& index, std::vector<size_t>& dim_list,
                               size_t dim_rem, size_t dim_sweep) {
    if (dim_rem == 0) {
      functor(source, result, index, dim_sweep);
    } else {
      typedef index_t index_type;

      AnovaHelper::level_t current_level;
      index_type current_index;

      index.get(dim_list[dim_rem - 1], current_level, current_index);

      if (current_level >= 1) {
        // given current point to next dim
        sweep_AnovaBoundary_rec(source, result, index, dim_list, dim_rem - 1, dim_sweep);

        if (!index.hint()) {
          index.leftChild(dim_list[dim_rem - 1]);

          if (!storage.isInvalidSequenceNumber(index.seq())) {
            sweep_AnovaBoundary_rec(source, result, index, dim_list, dim_rem, dim_sweep);
          }

          index.stepRight(dim_list[dim_rem - 1]);

          if (!storage.isInvalidSequenceNumber(index.seq())) {
            sweep_AnovaBoundary_rec(source, result, index, dim_list, dim_rem, dim_sweep);
          }

          index.up(dim_list[dim_rem - 1]);
        }
      } else {  
        // handle level minus one
        if (current_level == -1) {
          sweep_AnovaBoundary_rec(source, result, index, dim_list, dim_rem - 1, dim_sweep);
        }

        index.resetToLevelZeroInDim(dim_list[dim_rem - 1]);
        //Check if level zero exists
        if (!storage.isInvalidSequenceNumber(index.seq())) {
          sweep_AnovaBoundary_rec(source, result, index, dim_list, dim_rem - 1, dim_sweep);

          index.resetToLevelOneInDim(dim_list[dim_rem - 1]);

          // Check if level one exists
          if (!storage.isInvalidSequenceNumber(index.seq())) {
            sweep_AnovaBoundary_rec(source, result, index, dim_list, dim_rem, dim_sweep);
          }
        }

        index.resetToLevelMinusOneInDim(dim_list[dim_rem - 1]);
      }
    }
  }
};

}  // namespace base
}  // namespace sgpp

#endif /* SWEEP_HPP */
