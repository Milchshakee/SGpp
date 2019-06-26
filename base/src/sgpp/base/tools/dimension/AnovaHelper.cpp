#include "AnovaHelper.hpp"

namespace sgpp {
namespace base {
namespace AnovaHelper {

double getL2NormOfBasis(const AnovaHelper::LevelVector& levels) {
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

AnovaHelper::LevelVector getLevelVectorOfPoint(const GridPoint& point) {
  AnovaHelper::LevelVector v(point.getDimension());
  for (size_t d = 0; d < point.getDimension(); d++) {
    v[d] = point.getLevel(d);
  }
  return std::move(v);
}

  AnovaComponent getAnovaComponentOfPoint(const GridPoint& point) {
  AnovaComponent currentComp(point.getDimension(), false);
    for (size_t d = 0; d < point.getDimension(); d++) {
    currentComp[d] = point.getLevel(d) > 0;
  }
  
}

double calculateIncrementsVariance(GridStorage& gridStorage,
    const DataVector& alpha,
    const std::vector<LevelVector>& levels) {
  double sum = 0;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    GridPoint& gp = gridStorage.getPoint(i);
    LevelVector v = getLevelVectorOfPoint(gp);
    if (std::find(levels.begin(), levels.end(), v) != levels.end()) {
      double integral = getL2NormOfBasis(v);
      double val = alpha[i] * integral;
      sum += val;
    }
  }

  return sum;
}

double calculateIncrementVariance(GridStorage& gridStorage, const DataVector& alpha,
                                                              const AnovaHelper::LevelVector& level) {
  return calculateIncrementsVariance(gridStorage, alpha, {level});
}

bool isValidLevel(size_t maxLevel, const LevelVector& levels) {
  size_t sum = 0;
  for (size_t i = 0; i < levels.size(); i++) {
    sum += levels[i];
  }
  return sum <= maxLevel + levels.size() - 1;
}

bool nextLevelRaw(size_t maxLevel, LevelVector& levels,
                  const AnovaComponent& comp) {
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

bool nextLevel(size_t maxLevel, LevelVector& levels,
               const AnovaComponent& comp) {
  while (true) {
    bool b = nextLevelRaw(maxLevel, levels, comp);
    if (!b) {
      return false;
    }
    if (isValidLevel(maxLevel, levels)) {
      return true;
    }
  }
}

bool nextAnovaComponent(AnovaComponent& comp) {
  for (size_t i = 0; i < comp.size(); i++) {
    if (comp[i]) {
      comp[i] = false;
    } else {
      comp[i] = true;
      return true;
    }
  }
  return false;
}

double calculateDimensionVariance(GridStorage& gridStorage, const DataVector& alpha,
                                                              const AnovaHelper::AnovaComponent& comp) {
  AnovaHelper::LevelVector levels;
  for (size_t i = 0; i < comp.size(); i++) {
    if (!comp[i]) {
      levels.push_back(1);
    } else {
      levels.push_back(0);
    }
  }

  double sum = 0;
  while (nextLevel(gridStorage.getMaxLevel(), levels, comp)) {
    sum += calculateIncrementVariance(gridStorage, alpha, levels);
  }
  return sum;
}

Sample<AnovaComponent, double> calculateAnovaOrderVariances(GridStorage& gridStorage,
                                                            const DataVector& alpha) {
  AnovaComponent currentComp(gridStorage.getDimension(), false);
  std::vector<AnovaComponent> components;
  std::vector<double> variances;
  do {
    double var = calculateDimensionVariance(gridStorage, alpha, currentComp);
    components.push_back(currentComp);
    variances.push_back(var);
  } while (nextAnovaComponent(currentComp));
  return Sample<AnovaComponent, double>(components, variances);
}

  }
  }  // namespace base
}  // namespace sgpp
