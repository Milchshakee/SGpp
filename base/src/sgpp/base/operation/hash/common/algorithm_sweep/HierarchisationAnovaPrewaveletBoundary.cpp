// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/tools/sle/system/SLE.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/solver/Eigen.hpp>

namespace sgpp {
namespace base {

double BORDER_CENTER = 2.0;
double BORDER_OFFSET = 0.7;
double NORMAL_CENTER = 1.6;
double NORMAL_OFFSET = 0.4;

void HierarchisationAnovaPrewaveletBoundary::operator()(DataVector& source, DataVector& result,
                                                        grid_iterator& index, size_t dim) {
  AnovaBoundaryGrid::level_t max_level = index.getGridDepth(dim);

  if (max_level == -1) {
    return;
    }

  AnovaBoundaryGrid::level_t init_level;
  index_t init_index;
  index.get(dim, init_level, init_index);

  DataVector temp((1 << (max_level + 1)) + 2, 0);

  std::vector<index_t> lastIndices(max_level + 2, 0);
  for (AnovaBoundaryGrid::level_t l = 1; l <= max_level + 1; l++) {
    lastIndices.at(l) = (1 << l) - 1;
  }

  for (AnovaBoundaryGrid::level_t level = max_level; level >= 2; --level) {
    size_t n = (1 << (level - 1));
    DataMatrix m(n, n);
    m(0, 0) = BORDER_CENTER;
    m(0, 1) = BORDER_OFFSET;
    for (size_t i = 1; i <= n - 2; i++) {
      m(i, i) = NORMAL_CENTER;
      m(i, i - 1) = NORMAL_OFFSET;
      m(i, i + 1) = NORMAL_OFFSET;
    }
    m(n - 1, n - 1) = BORDER_CENTER;
    m(n - 1, n - 2) = BORDER_OFFSET;
    FullSLE sle(m);

    DataVector r(n, 0);
    DataVector u(n, 0);


    index_t firstIndex = 1;

    for (uint32_t i = 0; i < n; i++) {
      index_t ind = firstIndex + 2 * i;
      index.set(dim, level, ind);
      index_t nextLevelInd = 2 * ind;
      r.at(i) = result.at(index.seq()) - temp.at(nextLevelInd) + 0.5 * temp.at(nextLevelInd - 2) +
             0.5 * temp.at(nextLevelInd + 2);
    }

    sle_solver::Eigen a = sle_solver::Eigen();
    a.solve(sle, r, u);

    // Now the results are all ready in u.

    for (uint32_t i = 0; i < n; i++) {
      index.set(dim, level, (i * 2 + 1));
      result.at(index.seq()) = u.at(i);
    }

    DataVector newTemp((2 * n) + 2, 0);

    index.set(dim, level, 1);
    newTemp.at(0) = -1.2 * result.at(index.seq());
    if (level != max_level) {
      newTemp.at(0) += temp.at(0);
    }

    index.set(dim, level, lastIndices[level]);
    newTemp.at(lastIndices[level] + 1) = -1.2 * result.at(index.seq());
    if (level != max_level) {
      newTemp.at(lastIndices[level] + 1) += temp.at(lastIndices[level + 1] + 1);
    }

    // create new temp valuesand the normal start.
    // Please note, this is done in that strange
    for (index_t ind = 2; ind <= lastIndices.at(level) - 1; ind = ind + 2) {
      if (level != max_level) {
        newTemp.at(ind) = temp.at(ind * 2);
      }

      index.set(dim, level, ind - 1);
      double val = result.at(index.seq());
      newTemp.at(ind) -= 0.6 * val;

      index.set(dim, level, ind + 1);
      val = result.at(index.seq());
      newTemp.at(ind) -= 0.6 * val;
    }

    temp = newTemp;
  }

  double temp0 = 0.0;
  double temp2 = 0.0;
  if (max_level >= 1) {
    // Treatment of the top-point in this dimension
    index.resetToLevelOneInDim(dim);
    result.at(index.seq()) = 3.0 / 2.0 *
        (result.at(index.seq()) - temp.at(2) + (0.5 * temp.at(0)) + (0.5 * temp.at(4)));

    temp0 = result.at(index.seq()) + temp.at(0);
    temp2 = result.at(index.seq()) + temp.at(4);
  }

  index.resetToLevelMinusOneInDim(dim);
  double minusOneValue = result.at(index.seq());
  index.resetToLevelZeroInDim(dim);
  double zeroValue = result.at(index.seq());

  index.resetToLevelZeroInDim(dim);
  result.at(index.seq()) = (zeroValue - temp2 + temp0) / 2.0;
  index.resetToLevelMinusOneInDim(dim);
  result.at(index.seq()) = (2 * minusOneValue + zeroValue - temp2 - temp0) / 2.0;

  index.set(dim, init_level, init_index);
}

}  // namespace base
}  // namespace sgpp
