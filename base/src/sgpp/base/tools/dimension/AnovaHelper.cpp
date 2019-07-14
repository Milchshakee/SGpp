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

double getL2NormOfBasis(const GridPoint& gp) {
  size_t levelSum = 0;
  size_t level1Count = 0;
  for (size_t d = 0; d < gp.getDimension(); d++) {
    levelSum += std::max<size_t>(gp.getLevel(d), 1) - 1;
    if (gp.getLevel(d) == 1) {
      level1Count++;
      }
  }

  double f = std::pow(2.0 / 3.0, static_cast<double>(gp.getDimension()));
  double exp = std::pow(2.0, -(static_cast<double>(levelSum) - static_cast<double>(level1Count)));
  double val = std::sqrt(f * exp);

  return val;
}

AnovaComponent getAnovaComponentOfPoint(const GridPoint& point) {
  AnovaComponent currentComp(point.getDimension(), false);
  for (size_t d = 0; d < point.getDimension(); d++) {
    currentComp[d] = point.getLevel(d) > 0;
  }
  return currentComp;
}

Sample<AnovaComponent, double> calculateAnovaComponentVariances(GridStorage& gridStorage,
                                                            const DataVector& alpha) {
  std::map<AnovaComponent, double> variances;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    GridPoint& gp = gridStorage.getPoint(i);
    AnovaComponent c = getAnovaComponentOfPoint(gp);
    double integral = getL2NormOfBasis(gp);
    double val = std::abs(alpha[i]) * integral;
    variances.emplace(c, variances.find(c) == variances.end() ? val : variances.at(c) + val);
  }
  return Sample<AnovaComponent, double>(variances);
}

}  // namespace AnovaHelper
}  // namespace base
}  // namespace sgpp
