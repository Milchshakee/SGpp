#ifndef ANOVAHELPER
#define ANOVAHELPER

#include <vector>
#include "sgpp/base/grid/storage/hashmap/HashGridPoint.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {
namespace AnovaHelper {

class AnovaGridIterator {
 public:
  /// type of grid points
  typedef HashGridPoint index_type;
  /// index type
  typedef HashGridPoint::index_type index_t;
  /// level type
  typedef HashGridPoint::level_type level_t;

  explicit AnovaGridIterator(HashGridStorage& storage);

  /**
   *  Sets 0,0 in every dimension (Left Level zero ansatzfunction)
   */
  void resetToLevelZero();

  /**
   * left level zero ansatz function for a given dimension
   *
   * @param dim dimension in which we should step to level zero
   */
  void resetToLevelZeroInDim(size_t dim);

  /**
   * resets the iterator to the top if dimension d
   *
   * @param d the moving direction
   */
  void resetToLevelOneInDim(size_t d);

  /**
   * resets the iterator to the top if dimension d
   *
   * @param d the moving direction
   */
  void resetToLevelTwoInDim(size_t d);

  /**
   * left child in direction dim
   *
   * @param dim dimension in which we should step to the left child
   */
  void leftChild(size_t dim);

  /**
   * right child in direction dim
   *
   * @param dim dimension in which we should step to the right child
   */
  void rightChild(size_t dim);

  /**
   * hierarchical parent in direction dim
   *
   * @param d the moving direction
   */
  void up(size_t d);

  /**
   * step left in direction dim
   *
   * @param d the moving direction
   */
  void stepLeft(size_t d);

  /**
   * step right in direction dim
   *
   * @param d the moving direction
   */
  void stepRight(size_t d);

  /**
   * determines if the grid point is an inner grid point
   *
   * @return true if the grid point is an inner grid point
   */
  bool isInnerPoint() const;

  /**
   * returns true if there are no more children in any dimension
   *
   * @return returns true if there are no more children in any dimension
   */
  bool hint() const;

  /**
   * returns true if there are more left children in dimension d
   *
   * @param d the moving direction
   * @return true if there are more left children in dimension d
   */
  bool hintLeft(size_t d);

  /**
   * returns true if there are more right children in dimension d
   *
   * @param d the moving direction
   * @return true if there are more right children in dimension d
   */
  bool hintRight(size_t d);

  /**
   * Gets level @c l and index @c i in dimension @c d of the current grid point
   *
   * @param d the dimension of interest
   * @param l the ansatz function's level
   * @param i the ansatz function's index
   */
  inline void get(size_t d, index_type::level_type& l, index_type::index_type& i) const {
    index.get(d, l, i);
  }

  /**
   * Sets level @c l and index @c i in dimension @c d of the current grid point.
   * Recomputes the hash value of the current grid point.
   *
   * @param d the dimension of interest
   * @param l the ansatz function's level
   * @param i the ansatz function's index
   */
  inline void set(size_t d, index_type::level_type l, index_type::index_type i) {
    index.set(d, l, i);
    this->seq_ = storage.getSequenceNumber(index);
  }

  /**
   * Sets the iterator to a grid point.
   * Recomputes the hash value of the current grid point.
   *
   * @param point new grid point
   */
  inline void set(const index_type& point) {
    index = point;
    this->seq_ = storage.getSequenceNumber(index);
  }

  /**
   * Sets level @c l and index @c i in dimension @c d of the current grid point.
   * Does not recompute hash value of the current grid point.
   *
   * @param d the dimension of the gridpoint
   * @param l the ansatz function's level
   * @param i the ansatz function's index
   */
  inline void push(size_t d, index_type::level_type l, index_type::index_type i) {
    index.push(d, l, i);
  }

  /**
   * returns the current sequence number
   *
   * @return the current sequence number
   */
  size_t seq() const;

 private:
  /// reference the the hashmap that stores the gridpoints
  HashGridStorage& storage;
  /// GridPoint object used to operate on the current position in the hashmap
  HashGridPoint index;
  /// the current gridpoint's index
  size_t seq_;
};

/**
 * Vector that holds levels for every dimension
 */
typedef std::vector<sgpp::base::HashGridPoint::level_type> LevelVector;

typedef std::vector<bool> AnovaComponent;
typedef std::vector<AnovaComponent> AnovaComponentVector;

AnovaComponent getAnovaComponentOfPoint(const GridPoint& point);

Sample<AnovaComponent, double> calculateAnovaOrderVariances(GridStorage& gridStorage,
                                                            const DataVector& alpha);

}  // namespace AnovaHelper
}  // namespace base
}  // namespace sgpp

#endif