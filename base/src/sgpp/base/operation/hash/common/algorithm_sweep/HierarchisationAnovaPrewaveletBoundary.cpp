// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/tools/sle/system/SLE.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>

namespace sgpp {
namespace base {

void HierarchisationAnovaPrewaveletBoundary::hierarchiseBorders(
    DataVector& result, DataVector& functionValues, double& sum, size_t& count,
    AnovaGridIterator& originalIndex, size_t dim, AnovaGridIterator& index) {
  if (storage.isInvalidSequenceNumber(index.seq())) {
    return;
  }

  if (dim == storage.getDimension()) {
    double mult = 1;
    for (size_t d = 0; d < storage.getDimension(); d++) {
      index_t temp;
      AnovaTypes::level_t originalL;
      AnovaTypes::level_t l;
      originalIndex.get(d, originalL, temp);
      index.get(d, l, temp);
      if (originalL == 1 && l != -1) {
        mult *= -1;
      }

      if (originalL == 0 && l == -1) {
        mult *= -1;
      }
    }
    sum += functionValues.at(index.seq()) * mult;
    count++;
    return;
  }
  index.resetToLevelMinusOneInDim(dim);
  hierarchiseBorders(result, functionValues, sum, count, originalIndex, dim + 1, index);
  index.resetToLevelZeroInDim(dim);
  hierarchiseBorders(result, functionValues, sum, count, originalIndex, dim + 1, index);
  index.resetToLevelOneInDim(dim);
  hierarchiseBorders(result, functionValues, sum, count, originalIndex, dim + 1, index);
  index.resetToLevelMinusOneInDim(dim);
}

void HierarchisationAnovaPrewaveletBoundary::convertLevelZeroAnsatzFunction(DataVector& source, DataVector& result, AnovaGridIterator& index, size_t dim) {
  if (dim == storage.getDimension()) {
    AnovaGridIterator it(storage);
    double sum = 0;
    size_t count = 0;
    hierarchiseBorders(result, source, sum, count, index, 0, it);
    result.at(index.seq()) = sum / count;
    return;
    }

  
  index.resetToLevelMinusOneInDim(dim);
    convertLevelZeroAnsatzFunction(source, result, index, dim + 1);
  index.resetToLevelZeroInDim(dim);
  convertLevelZeroAnsatzFunction(source, result, index, dim + 1);
  index.resetToLevelOneInDim(dim);
  convertLevelZeroAnsatzFunction(source, result, index, dim + 1);
  index.resetToLevelMinusOneInDim(dim);
}

void HierarchisationAnovaPrewaveletBoundary::convertLevelZeroAnsatzFunction(DataVector& source,
                                                                            DataVector& result) {
  AnovaGridIterator it(storage);
  convertLevelZeroAnsatzFunction(source, result, it, 0);
}

double BORDER_CENTER = 2.0;
double BORDER_OFFSET = 0.7;
double NORMAL_CENTER = 1.6;
double NORMAL_OFFSET = 0.4;

void HierarchisationAnovaPrewaveletBoundary::operator()(DataVector& source, DataVector& result,
                                                        grid_iterator& index, size_t dim) {
  AnovaBoundaryGrid::level_t max_level = index.getGridDepth(dim);

  if (max_level <= 1) {
    return;
  }

  AnovaBoundaryGrid::level_t init_level;
  index_t init_index;
  index.get(dim, init_level, init_index);

  std::vector<double> temp(1 << (max_level + 1), 0);
  std::vector<double> r(1 << (max_level - 1), 0);
  std::vector<double> gam(1 << (max_level - 1), 0);
  DataVector u(1 << (max_level - 1), 0);

  std::vector<index_t> lastIndices(max_level + 2, 0);
  for (AnovaBoundaryGrid::level_t l = 1; l <= max_level + 1; l++) {
    lastIndices.at(l) = (1 << l) - 1;
  }

  for (AnovaBoundaryGrid::level_t level = max_level; level >= 2; --level) {
    size_t n = (1 << (level));
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



    index_t firstIndex = 1;
    size_t equationSystemSize = (1 << (level - 1));

    // First, we set the right-hand side of the triangular eqaution system

    // special treatment for the first point in the row
    index.set(dim, level, firstIndex);
    r.at(0) = result.at(index.seq()) - temp.at(firstIndex) + 0.5 * temp.at(firstIndex + 2);

    // special treatment for the last point in the row
    index.set(dim, level, lastIndices.at(level));
    r.at(equationSystemSize - 1) = result.at(index.seq()) - temp.at(lastIndices.at(level + 1) - 2) +
                                0.5 * temp.at(lastIndices.at(level + 1) - 4);

    // normal treatment for the middle
    if (level > 2) {
      for (uint32_t i = 1; i < equationSystemSize - 1; i++) {
        index_t ind = firstIndex + 2 * i;
        index.set(dim, level, ind);
        index_t nextLevelInd = firstIndex + 4 * i;
        r.at(i) = result.at(index.seq()) - temp.at(nextLevelInd) + 0.5 * temp.at(nextLevelInd - 2) +
               0.5 * temp.at(nextLevelInd + 2);
      }
    }

        // Run the actual triangulation

    double bet = 0.0;

    // This is the forward-reduction
    for (int i = 0; i < 1 << (level - 1); i++) {
      index.set(dim, level, i * 2 + 1);

      if (storage.isInvalidSequenceNumber(index.seq())) {
        u[i] = 0;
        gam[i] = 0;
        continue;
      } else {
        if (i == 0) {
          bet = 1.2;
          u[0] = (r[0] / bet);
          continue;
        }

        /* We need to distinguish between the start of the triangulation
         * and the normal start. Please note, this is done in that strange
         * because in the adaptive case it is perfectly possible to have
         * grid points just in the middle of the level, thus the beginning
         * of the triangulation is not determined with index 1.
         */

        index.set(dim, level, i * 2 - 1);

        if (!storage.isInvalidSequenceNumber(index.seq()))
          gam[i] = 0.4 / bet;
        else
          gam[i] = 0.0;

        if (i == (1 << (level - 1)) - 1) {
          bet = 1.2;
        } else {
          bet = 1.6;
        }

        index.set(dim, level, i * 2 + 1);
        bet = bet - 0.4 * gam[i];
        u[i] = (r[i] - 0.4 * u[i - 1]) / bet;
      }
    }

    //// Run the actual triangulation

    //u.at(0) = (r.at(0) / 1.2);

    //double lastT = 1.6;
    //// This is the forward-reduction
    //for (index_t i = 1; i < lastIndices.at(level - 1); i++) {
    //  gam.at(i) = 0.4 / lastT;

    //  lastT = 1.6 - 0.4 * gam.at(i);
    //  u.at(i) = (r.at(i) - 0.4 * u.at(i - 1)) / lastT;
    //}

    //index_t l = lastIndices.at(level - 1);
    //gam.at(l) = 0.4 / 1.2;
    //u.at(l) = (r.at(l) - 0.4 * u.at(l - 1)) / (1.2 - 0.4 * gam.at(l));

    // Backward-Reduction
    for (int i = lastIndices.at(level - 1) - 1; i >= 0; --i) {
      u.at(i) = u.at(i) - gam.at(i + 1) * u.at(i + 1);
    }

    // Now the results are all ready in u.

    for (uint32_t i = 0; i <= lastIndices.at(level - 1); i++) {
      index.set(dim, level, (i * 2 + 1));
      result.at(index.seq()) = u.at(i);
    }

    // create new temp valuesand the normal start.
    // Please note, this is done in that strange
    for (index_t ind = 1; ind <= lastIndices.at(level) - 2; ind = ind + 2) {
      index.set(dim, level, ind);
      double val = result.at(index.seq());

      if (level != max_level) {
        temp.at(ind) = temp.at(ind * 2 + 1);
      }

      temp.at(ind) = temp.at(ind) - 0.6 * val;

      index.set(dim, level, ind + 2);
      val = result.at(index.seq());

      temp.at(ind) = temp.at(ind) - 0.6 * val;
    }
  }

  // Treatment of the top-point in this dimension
  index.resetToLevelOneInDim(dim);
  result.at(index.seq()) -= temp.at(1);

  index.set(dim, init_level, init_index);
}

}  // namespace base
}  // namespace sgpp
