#include "AnovaHelper.hpp"

namespace sgpp {
namespace base {
namespace AnovaHelper {

AnovaGridIterator::AnovaGridIterator(HashGridStorage& storage)
    : storage(storage), index(storage.getDimension()) {
  resetToLevelZero();
  }

void AnovaGridIterator::resetToLevelZero() {
  for (size_t i = 0; i < storage.getDimension(); i++) {
    index.push(i, 0, 0);
  }

  index.rehash();
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::resetToLevelZeroInDim(size_t dim) {
  index.set(dim, 0, 0);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::resetToLevelOneInDim(size_t d) {
  index.set(d, 1, 2);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::resetToLevelTwoInDim(size_t d) {
  index.set(d, 2, 2);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::leftChild(size_t dim) {
  HashGridPoint::level_type l;
  HashGridPoint::index_type i;
  index.get(dim, l, i);
  index.set(dim, l + 1, 2 * (i - 1));
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::rightChild(size_t dim) {
  HashGridPoint::level_type l;
  HashGridPoint::index_type i;
  index.get(dim, l, i);
  index.set(dim, l + 1, 2 * (i + 1));
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::up(size_t d) {
  HashGridPoint::level_type l;
  HashGridPoint::index_type i;
  index.get(d, l, i);
  i /= 2;
  i = i % 4 == 1 ? i + 1 : i - 1;
  index.set(d, l - 1, i);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::stepLeft(size_t d) {
  HashGridPoint::level_type l;
  HashGridPoint::index_type i;
  index.get(d, l, i);
  index.set(d, l, i - 4);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::stepRight(size_t d) {
  HashGridPoint::level_type l;
  HashGridPoint::index_type i;
  index.get(d, l, i);
  index.set(d, l, i + 4);
  this->seq_ = storage.getSequenceNumber(index);
}

bool AnovaGridIterator::isInnerPoint() const { return index.isInnerPoint(); }

bool AnovaGridIterator::hint() const { return storage.getPoint(this->seq_).isLeaf(); }

bool AnovaGridIterator::hintLeft(size_t d) {
  HashGridPoint::level_type l;
  HashGridPoint::index_type i;
  index.get(d, l, i);
  index.set(d, l + 1, 2 * (i - 1));

  bool hasIndex = storage.isContaining(index);

  index.set(d, l, i);

  return hasIndex;
}

bool AnovaGridIterator::hintRight(size_t d) {
  HashGridPoint::level_type l;
  HashGridPoint::index_type i;
  index.get(d, l, i);
  index.set(d, l + 1, 2 * (i + 1));

  bool hasIndex = storage.isContaining(index);

  index.set(d, l, i);

  return hasIndex;
}

size_t AnovaGridIterator::seq() const { return seq_; }

double getL2NormOfBasis(const AnovaHelper::LevelVector& levels) {
  size_t levelSum = 0;
  for (size_t d = 0; d < levels.size(); d++) {
    levelSum += levels[d] - 1;
  }
  return std::pow(2.0 / 3.0, static_cast<double>(levels.size()) / 2.0) / std::pow(2.0, static_cast<double>(levelSum) / 2.0);
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
  return currentComp;
}

double calculateIncrementsVariance(GridStorage& gridStorage, const DataVector& alpha,
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

bool nextLevelRaw(size_t maxLevel, LevelVector& levels, const AnovaComponent& comp) {
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

bool nextLevel(size_t maxLevel, LevelVector& levels, const AnovaComponent& comp) {
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

}  // namespace AnovaHelper
}  // namespace base
}  // namespace sgpp
