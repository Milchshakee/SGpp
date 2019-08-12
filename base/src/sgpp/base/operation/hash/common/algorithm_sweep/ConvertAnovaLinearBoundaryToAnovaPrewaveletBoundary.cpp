// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary.hpp>

namespace sgpp {
namespace base {

void ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary::convertLevelZeroAnsatzFunction(
  DataVector& source, DataVector& result, grid_iterator& index, size_t dim) {
  index.resetToLevelMinusOneInDim(dim);
  size_t constantSeq = index.seq();

  index.resetToLevelZeroInDim(dim);
  if (!storage.isInvalidSequenceNumber(index.seq())) {
    result[constantSeq] += source[index.seq()] / 2;
    result[index.seq()] -= source[index.seq()] / 2;
  }
}

void ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary::operator()(DataVector& source,
                                                               DataVector& result,
                                           grid_iterator& index, size_t dim) {
  convertLevelZeroAnsatzFunction(source, result, index, dim);
  index.resetToLevelOneInDim(dim);
  if (!storage.isInvalidSequenceNumber(index.seq())) {
    HashGridIterator it(storage);
    it.set(index.getIndex());
    convert.operator()(source, result, it, dim);
  }
}

}  // namespace base
}  // namespace sgpp
