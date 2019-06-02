#include <sgpp/base/operation/hash/OperationAnova.hpp>
#include "common/basis/AnovaBoundaryBasis.hpp"
#include "sgpp/base/grid/GridStorage.hpp"
#include "sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp"

double getL2NormOfBasis(const sgpp::base::OperationAnova::LevelVector& levels) {
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

sgpp::base::OperationAnova::LevelVector getLevelVector(const sgpp::base::GridPoint& point) {
  sgpp::base::OperationAnova::LevelVector v(point.getDimension());
  for (size_t d = 0; d < point.getDimension(); d++) {
    v[d] = point.getLevel(d);
  }
  return std::move(v);
}

double sgpp::base::OperationAnova::calculateIncrementsVariance(
    const sgpp::base::DataVector& alpha,
    const std::vector<sgpp::base::OperationAnova::LevelVector>& levels) {
  double sum = 0;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::OperationAnova::LevelVector v = getLevelVector(gp);
    if (std::find(levels.begin(), levels.end(), v) != levels.end()) {
      double integral = getL2NormOfBasis(v);
      double val = alpha[i] * integral;
      sum += val;
    }
  }

  return sum;
}

double sgpp::base::OperationAnova::calculateIncrementVariance(const DataVector& alpha,
                                                                 const LevelVector& level) {
  return calculateIncrementsVariance(alpha, {level});
}

bool isValidLevel(const sgpp::base::OperationAnova::LevelVector& levels) { return true; }

bool nextLevel(size_t maxLevel, sgpp::base::OperationAnova::LevelVector& levels,
               const sgpp::base::OperationAnova::DimensionVector& fixedDimensions) {
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

double sgpp::base::OperationAnova::calculateDimensionVariance(
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

bool nextDimensionVector(
    sgpp::base::OperationAnova::DimensionVector& fixedDimensions) {
  size_t first = 0;
  for (size_t i = 0; i < fixedDimensions.size(); i++) {
    if (!fixedDimensions[i]) {
      first = i;
      fixedDimensions[i] = true;
      break;
      }
  }

  size_t propagated = 0;
    for (size_t i = first + 1; i < fixedDimensions.size(); i++) {
    if (fixedDimensions[i]) {
      fixedDimensions[i] = false;
      break;
    } else {
      fixedDimensions[i] = true;
      propagated++;
    }
  }

  if (first + propagated + 1 == fixedDimensions.size()) {
    return false;
    }

    for (size_t i = 0; i < propagated; i++) {
    fixedDimensions[i] = false;
  }

  return true;
}


std::vector<sgpp::base::OperationAnova::AnovaComponent> sgpp::base::OperationAnova::
calculateAnovaOrderVariance(const DataVector& alpha, size_t anovaOrder) {
  DimensionVector fixedDimensions(gridStorage.getDimension(), false);
  for (size_t i = 0; i < anovaOrder; i++) {
    fixedDimensions[i] = true;
  }

  std::vector<sgpp::base::OperationAnova::AnovaComponent> components;
  while (nextDimensionVector(fixedDimensions))
  {
    AnovaComponent comp{anovaOrder, fixedDimensions,
                        calculateDimensionVariance(alpha, fixedDimensions)};
    components.push_back(comp);
  }
  return std::move(components);
}
