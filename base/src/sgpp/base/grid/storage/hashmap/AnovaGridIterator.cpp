
#include <sgpp/base/grid/storage/hashmap/AnovaGridIterator.hpp>

namespace sgpp {
namespace base {

AnovaGridIterator::AnovaGridIterator(HashGridStorage& storage)
    : storage(storage), index(storage.getDimension()) {
  resetToLevelMinusOne();
}

const HashGridPoint& AnovaGridIterator::getIndex() { return index; }

void AnovaGridIterator::resetToLevelMinusOne() {
  for (size_t i = 0; i < storage.getDimension(); i++) {
    index.push(i, 0, 0);
  }

  index.rehash();
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::resetToLevelMinusOneInDim(size_t dim) {
  index.set(dim, 0, 0);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::resetToLevelZeroInDim(size_t d) {
  index.set(d, 0, 1);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::resetToLevelOneInDim(size_t d) {
  index.set(d, 1, 1);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::leftChild(size_t dim) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(dim, l, i);
  index.set(dim, l + 1, 2 * i - 1);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::rightChild(size_t dim) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(dim, l, i);
  index.set(dim, l + 1, 2 * i + 1);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::up(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(d, l, i);

  i /= 2;
  i += i % 2 == 0 ? 1 : 0;

  index.set(d, l - 1, i);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::stepLeft(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(d, l, i);
  index.set(d, l, i - 2);
  this->seq_ = storage.getSequenceNumber(index);
}

void AnovaGridIterator::stepRight(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(d, l, i);
  index.set(d, l, i + 2);
  this->seq_ = storage.getSequenceNumber(index);
}

bool AnovaGridIterator::isInnerPoint() const { return index.isInnerPoint(); }

bool AnovaGridIterator::hint() const { return storage.getPoint(this->seq_).isLeaf(); }

bool AnovaGridIterator::hintLeft(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  bool hasIndex = true;

  index.get(d, l, i);
  index.set(d, l + 1, 2 * i - 1);

  hasIndex = storage.isContaining(index);

  index.set(d, l, i);

  return hasIndex;
}

bool AnovaGridIterator::hintRight(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  bool hasIndex = true;

  index.get(d, l, i);
  index.set(d, l + 1, 2 * i + 1);

  hasIndex = storage.isContaining(index);

  index.set(d, l, i);

  return hasIndex;
}

void AnovaGridIterator::get(size_t d, AnovaBoundaryGrid::level_t& l, index_t& i) const {
  HashGridPoint::level_type lRaw;
  HashGridPoint::index_type iRaw;
  index.get(d, lRaw, iRaw);
  AnovaBoundaryGrid::fromNormalGridPointLevelIndex(lRaw, iRaw, l, i);
}

void AnovaGridIterator::set(size_t d, AnovaBoundaryGrid::level_t l, index_t i) {
  HashGridPoint::level_type lRaw;
  HashGridPoint::index_type iRaw;
  AnovaBoundaryGrid::toNormalGridPointLevelIndex(l, i, lRaw, iRaw);
  index.set(d, lRaw, iRaw);
}

AnovaBoundaryGrid::level_t AnovaGridIterator::getGridDepth(size_t dim) {
  AnovaBoundaryGrid::level_t depth = -1;

  AnovaBoundaryGrid::level_t orig_level, cur_level;
  index_t orig_index, cur_index;
  get(dim, orig_level, orig_index);

  resetToLevelMinusOneInDim(dim);
  resetToLevelZeroInDim(dim);
  if (!storage.isInvalidSequenceNumber(this->seq())) {
    depth++;
    resetToLevelOneInDim(dim);
    if (!storage.isInvalidSequenceNumber(this->seq())) {
      depth++;

      while (true) {
        if (this->hintLeft(dim)) {
          depth++;
          this->leftChild(dim);
        } else if (this->hintRight(dim)) {
          depth++;
          this->rightChild(dim);
        } else {
          break;
        }
      }
    }
  }

  this->set(dim, orig_level, orig_index);
  return depth;
}

size_t AnovaGridIterator::seq() const { return seq_; }

}  // namespace base
}  // namespace sgpp
