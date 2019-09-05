// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary.hpp>

namespace sgpp {
namespace base {

double ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary::calcBorderFuntionValues(
    DataVector& source, DataVector& functionValues, AnovaBoundaryGrid::AnovaComponent& comp) {
  AnovaGridIterator index(storage);
  for (size_t d = 0; d < storage.getDimension(); d++) {
    if (comp[d]) {
      index.set(d, 0, 1);
      AnovaBoundaryGrid::AnovaComponent copy = comp;
      copy[d] = 0;
      functionValues[index.seq()] += calcBorderFuntionValues(source, functionValues, copy);
    } else {
      index.set(d, -1, 0);
    }
  }

  return functionValues[index.seq()];
}

double ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary::hierarchiseBorders(
    DataVector& functionValues, AnovaBoundaryGrid::AnovaComponent& comp, bool neg) {
  AnovaGridIterator index(storage);
  double res = 0;
  for (size_t d = 0; d < storage.getDimension(); d++) {
    if (comp[d]) {
      index.set(d, 0, 1);
      AnovaBoundaryGrid::AnovaComponent copy = comp;
      copy[d] = 0;
      res += hierarchiseBorders(functionValues, copy, !neg);
    } else {
      index.set(d, -1, 0);
    }
  }
  if (neg) {
    res -= functionValues[index.seq()];
  } else {
    res += functionValues[index.seq()];
  }

  return res / 2.0;
}

double ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary::hierarchiseBorders(
    DataVector& result, DataVector& functionValues) {
  AnovaBoundaryGrid::AnovaComponent comp(storage.getDimension());
  AnovaGridIterator index(storage);
  for (size_t i = 1; i <= 1 << storage.getDimension(); ++i) {
    for (size_t d = 0; d < storage.getDimension(); d++) {
      if (i & (1 << d)) {
        comp[i] = true;
        index.set(d, 0, 1);
      } else {
        comp[i] = false;
        index.set(d, -1, 0);
      }

      result[index.seq()] = hierarchiseBorders(functionValues, comp, false);
    }
  }
}

void ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary::convertLevelZeroAnsatzFunction(
    DataVector& source, DataVector& result) {
  AnovaGridIterator it(storage);
  DataVector funcValues(source.size(), 0);
  auto c = AnovaBoundaryGrid::AnovaComponent(storage.getDimension(), true);
  funcValues[0] = source[0];
  calcBorderFuntionValues(source, funcValues, c);
  hierarchiseBorders(result, funcValues);
}

void ConvertAnovaLinearBoundaryToAnovaPrewaveletBoundary::operator()(DataVector& source,
                                                                     DataVector& result,
                                                                     grid_iterator& index,
                                                                     size_t dim) {
  // The Algorithm runs bottom-up, thus we need the maximal grid-depth
  // in the current dimension
  AnovaBoundaryGrid::level_t max_level = index.getGridDepth(dim);

  if (max_level == -1) {
    return;
  }

  AnovaBoundaryGrid::level_t init_level;
  index_t init_index;
  index.get(dim, init_level, init_index);

  double leftBorder = source[index.seq()];
  index.resetToLevelZeroInDim(dim);
  double rightBorder = source[index.seq()];
  index.resetToLevelOneInDim(dim);
  double middle = storage.isInvalidSequenceNumber(index.seq()) ? 0.0 : source[index.seq()];

  index.resetToLevelMinusOneInDim(dim);
  result[index.seq()] = (leftBorder + rightBorder + 2 * middle) / 4.0;
  index.resetToLevelZeroInDim(dim);
  result[index.seq()] = (-leftBorder + rightBorder) / 2.0;
  index.resetToLevelOneInDim(dim);
  if (!storage.isInvalidSequenceNumber(index.seq())) {
    result[index.seq()] = (-leftBorder - rightBorder + 2 * middle) / 4.0;
  }

  if (max_level <= 1) {
    return;
  }

  std::vector<double> temp(1 << (max_level + 1), 0);
  std::vector<double> r(1 << (max_level - 1), 0);
  std::vector<double> gam(1 << (max_level - 1), 0);
  std::vector<double> u(1 << (max_level - 1), 0);

  std::vector<index_t> lastIndices(max_level + 1, 0);
  for (size_t l = 1; l <= max_level; l++) {
    lastIndices[l] = (1 << l) - 1;
  }

  for (AnovaBoundaryGrid::level_t level = max_level; level >= 2; --level) {
    index_t firstIndex = 1;
    size_t equationSystemSize = (1 << (level - 1));

    // First, we set the right-hand side of the triangular eqaution system

    // special treatment for the first point in the row
    index.set(dim, level, firstIndex);
    r[0] = result[index.seq()] - temp[firstIndex] + 0.5 * temp[firstIndex + 2];

    // special treatment for the last point in the row
    index.set(dim, level, lastIndices[level]);
    r[equationSystemSize - 1] = result[index.seq()] - temp[lastIndices[level + 1] - 2] +
                                0.5 * temp[lastIndices[level + 1] - 4];

    // normal treatment for the middle
    if (level > 2) {
      for (size_t i = 1; i < equationSystemSize - 1; i++) {
        index_t ind = firstIndex + 2 * i;
        index.set(dim, level, ind);
        index_t nextLevelInd = firstIndex + 4 * i;
        r[i] = result[index.seq()] - temp[nextLevelInd] + 0.5 * temp[nextLevelInd - 2] +
               0.5 * temp[nextLevelInd + 2];
      }
    }

    // Run the actual triangulation

    double bet = 1.2;
    u[0] = (r[0] / 1.2);

    // This is the forward-reduction
    for (int i = 1; i < lastIndices[level - 1]; i++) {
      index_t ind = firstIndex + i * 2;
      gam[i] = 0.4 / bet;

      bet = 1.6 - 0.4 * gam[i];
      u[i] = (r[i] - 0.4 * u[i - 1]) / bet;
    }

    size_t l = lastIndices[level - 1];
    gam[l] = 0.4 / bet;
    u[l] = (r[l] - 0.4 * u[l - 1]) / (1.2 - 0.4 * gam[l]);

    // Backward-Reduction
    for (int i = lastIndices[level - 1] - 1; i >= 0; --i) {
      u[i] = u[i] - gam[i + 1] * u[i + 1];
    }

    // Now the results are all ready in u.

    for (int i = 0; i <= lastIndices[level - 1]; i++) {
      index.set(dim, level, (i * 2 + 1));
      result[index.seq()] = u[i];
    }

    // create new temp valuesand the normal start.
    // Please note, this is done in that strange
    for (index_t ind = 1; ind <= lastIndices[level]; ind = ind + 2) {
      index.set(dim, level, ind);
      double val = result[index.seq()];

      if (level != max_level) {
        temp[ind] = temp[ind * 2 + 1];
      }

      temp[ind] = temp[ind] - 0.6 * val;

      index.set(dim, level, ind + 2);
      val = result[index.seq()];

      temp[ind] = temp[ind] - 0.6 * val;
    }
  }

  // Treatment of the top-point in this dimension
  index.set(dim, init_level, init_index);
  result[index.seq()] -= temp[1];
}

}  // namespace base
}  // namespace sgpp
