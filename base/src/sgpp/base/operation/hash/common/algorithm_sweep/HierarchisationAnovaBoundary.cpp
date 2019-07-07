// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

HierarchisationAnovaBoundary::HierarchisationAnovaBoundary(GridStorage& storage)
    : storage(storage) {}

HierarchisationAnovaBoundary::~HierarchisationAnovaBoundary() {}

void HierarchisationAnovaBoundary::hierarchiseConstantRec(DataVector& source, DataVector& result,
                                                          grid_iterator& index, size_t dim,
                                                          double constant) {
  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    // descend left
    index.leftChild(dim);
    if (!storage.isInvalidSequenceNumber(index.seq())) {
      hierarchiseConstantRec(source, result, index, dim, constant);
    }

    // descend right
    index.stepRight(dim);
    if (!storage.isInvalidSequenceNumber(index.seq())) {
      hierarchiseConstantRec(source, result, index, dim, constant);
    }

    // ascend
    index.up(dim);
  }
  // hierarchisation
  size_t seq = index.seq();
  result[seq] -= constant;
}

void HierarchisationAnovaBoundary::hierarchiseConstant(DataVector& source, DataVector& result,
                                                       grid_iterator& index, size_t dim) {
  // left constant boundary
  index.resetToLevelZeroInDim(dim);
  double constant = source[index.seq()];

  index.resetToLevelOneInDim(dim);
  if (!storage.isInvalidSequenceNumber(index.seq())) {
    hierarchiseConstantRec(source, result, index, dim, constant);
  }
}



void HierarchisationAnovaBoundary::rec(DataVector& source, DataVector& result, grid_iterator& index,
                                       size_t dim, double fl, double fr) {
  // current position on the grid
  size_t seq = index.seq();
  // value in the middle, needed for recursive call and
  // calculation of the hierarchical surplus
  double fm = source[seq];

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    // descend left
    index.leftChild(dim);
    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fl, fm);
    }

    // descend right
    index.stepRight(dim);
    if (!storage.isInvalidSequenceNumber(index.seq())) {
      rec(source, result, index, dim, fm, fr);
    }

    // ascend
    index.up(dim);
  }

  // hierarchisation
  result[seq] = fm - ((fl + fr) / 2.0);
}

void HierarchisationAnovaBoundary::operator()(DataVector& source, DataVector& result,
                                              grid_iterator& index, size_t dim) {
  hierarchiseConstant(source, result, index, dim);

  index.resetToLevelOneInDim(dim);
  if (!storage.isInvalidSequenceNumber(index.seq())) {
    double right_boundary = source[index.seq()];

      // see if we can go down further
    if (!index.hint()) {
      index.resetToLevelTwoInDim(dim);
      if (!storage.isInvalidSequenceNumber(index.seq())) {
        rec(source, result, index, dim, 0, right_boundary);
      }
    }
  }
  index.resetToLevelZeroInDim(dim);
}

}  // namespace base
}  // namespace sgpp
