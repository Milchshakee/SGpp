#include <sgpp/base/operation/hash/OperationVariance.hpp>
#include "common/basis/AnovaBoundaryBasis.hpp"
#include "sgpp/base/grid/GridStorage.hpp"
#include "sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp"

double getL2NormOfBasis(const sgpp::base::OperationVariance::LevelVector& levels) {
  double result = 1;
  for (size_t d = 0; d < levels.size(); d++) {
    if (levels[d] == 0) {
      result *= 1;
    } else {
      double integral = (1. / 3.) / static_cast<double>(1 << (levels[d] - 2));
      result *= integral;
    }
  }
  return result;
}

sgpp::base::OperationVariance::LevelVector getLevelVector(const sgpp::base::GridPoint& point) {
  sgpp::base::OperationVariance::LevelVector v(point.getDimension());
  for (size_t d = 0; d < point.getDimension(); d++) {
    v[d] = point.getLevel(d);
  }
  return std::move(v);
}

double sgpp::base::OperationVariance::calculateIncrementsVariance(
    const sgpp::base::DataVector& alpha,
    const std::vector<sgpp::base::OperationVariance::LevelVector>& levels) {
  double sum = 0;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::OperationVariance::LevelVector v = getLevelVector(gp);
    if (std::find(levels.begin(), levels.end(), v) != levels.end()) {
      double integral = getL2NormOfBasis(v);
      double val = alpha[i] * integral;
      sum += val;
    }
  }

  return sum;
}

double sgpp::base::OperationVariance::calculateIncrementVariance(const DataVector& alpha,
                                                                 const LevelVector& level) {
  return calculateIncrementsVariance(alpha, {level});
}

bool isValidLevel(const sgpp::base::OperationVariance::LevelVector& levels) { return true; }

bool nextLevel(size_t maxLevel, sgpp::base::OperationVariance::LevelVector& levels,
               const sgpp::base::OperationVariance::DimensionVector& fixedDimensions) {
  bool lastAdded = false;
  for (size_t i = 0; i < fixedDimensions.size(); i++) {
    if (!fixedDimensions[i]) {
      if (levels[i] == maxLevel) {
        levels[i] = 0;
        lastAdded = true;
      } else {
        lastAdded = false;
        levels[i]++;
      }

      if (!lastAdded) {
        return true;
      }
    }
  }

  // Reached end
  return false;
}

double sgpp::base::OperationVariance::calculateDimensionVariance(
    const DataVector& alpha, const DimensionVector& fixedDimensions) {
  LevelVector dimensions;
  for (size_t i = 0; i < fixedDimensions.size(); i++) {
    if (!fixedDimensions[i]) {
      dimensions.push_back(1);
    } else {
      dimensions.push_back(0);
    }
  }

  double sum = 0;
  while (nextLevel(gridStorage.getMaxLevel(), dimensions, fixedDimensions)) {
    if (isValidLevel(dimensions)) {
      sum += calculateIncrementVariance(alpha, dimensions);
    }
  }
  return sum;
}
