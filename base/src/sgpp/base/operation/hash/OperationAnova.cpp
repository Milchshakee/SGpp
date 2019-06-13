#include <sgpp/base/operation/hash/OperationAnova.hpp>
#include "common/basis/AnovaBoundaryBasis.hpp"
#include "sgpp/base/grid/GridStorage.hpp"
#include "sgpp/base/tools/Sample.hpp"

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
               const sgpp::base::AnovaComponent& comp) {
  bool lastAdded = false;
  for (size_t i = 0; i < comp.size(); i++) {
    if (!comp[i]) {
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

bool nextAnovaComponent(
    sgpp::base::AnovaComponent& comp) {
  size_t first = 0;
  for (size_t i = 0; i < comp.size(); i++) {
    if (!comp[i]) {
      first = i;
      comp[i] = true;
      break;
      }
  }

  size_t propagated = 0;
    for (size_t i = first + 1; i < comp.size(); i++) {
    if (comp[i]) {
      comp[i] = false;
      break;
    } else {
      comp[i] = true;
      propagated++;
    }
  }

  if (first + propagated + 1 == comp.size()) {
    return false;
    }

    for (size_t i = 0; i < propagated; i++) {
    comp[i] = false;
  }

  return true;
}

double sgpp::base::OperationAnova::calculateDimensionVariance(const DataVector& alpha,
                                                              const AnovaComponent& comp) {
  LevelVector levels;
  for (size_t i = 0; i < comp.size(); i++) {
    if (!comp[i]) {
      levels.push_back(1);
    } else {
      levels.push_back(0);
    }
  }

  double sum = 0;
  while (nextLevel(gridStorage.getMaxLevel(), levels, comp)) {
    if (isValidLevel(levels)) {
      sum += calculateIncrementVariance(alpha, levels);
    }
  }
  return sum;
}

sgpp::base::Sample<sgpp::base::AnovaComponent, double> sgpp::base::OperationAnova::calculateAnovaOrderVariances(
    const DataVector& alpha) {
  AnovaComponent currentComp(gridStorage.getDimension(), false);
  std::vector<AnovaComponent> components;
  std::vector<double> variances;
  while (nextAnovaComponent(currentComp)) {
    double var = calculateDimensionVariance(alpha, currentComp);
    components.push_back(currentComp);
    variances.push_back(var);
  }
  return Sample<AnovaComponent, double>(components, variances);
}
