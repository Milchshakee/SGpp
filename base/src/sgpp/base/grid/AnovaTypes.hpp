#pragma once

#include <cstdint>
#include <sgpp/base/grid/storage/hashmap/HashGridPoint.hpp>

namespace sgpp {
namespace base {
namespace AnovaTypes {

/**
 * Level type for ANOVA grids (uses -1 level)
 */
typedef int32_t level_t;

  /**
 * Level type for ANOVA grids (uses -1 level)
 */
typedef uint32_t index_t;

/**
 * Level-Index pair for ANOVA grids (uses -1 level)
 */
struct LevelIndexPair {
  /**
   * Level of the grid point in the hierarchy.
   */
  level_t level;
  /**
   * Index of the grid point in the index set for that particular level.
   */
  sgpp::base::HashGridPoint::index_type index;
};

    /**
 * Converts a compatibility level to the true ANOVA grid level
 */
inline level_t fromNormalLevel(base::level_t l) {
  if (l == 0) {
    return -1;
  } else {
    return static_cast<level_t>(l - 1);
  }
}

/**
 * Converts an ANOVA grid level to the compatibility level
 */
inline base::level_t toNormalLevel(level_t l) {
  if (l == -1) {
    return 0;
  } else {
    return static_cast<base::level_t>(l + 1);
  }
}

/**
 * Converts the level and index ofa true ANOVA grid point to the level and index of a compatibility
 * grid point
 */
inline void toNormalGridPointLevelIndex(level_t l, HashGridPoint::index_type i,
                                        HashGridPoint::level_type& lOut,
                                        HashGridPoint::index_type& iOut) {
  if (l == -1) {
    lOut = 0;
    iOut = 0;
  } else if (l == 0) {
    lOut = 0;
    iOut = 1;
  } else {
    lOut = l;
    iOut = i;
  }
}

/**
 * Converts the level and index of a compatibility grid point to the level and index of a true
 * ANOVA grid point
 */
inline void fromNormalGridPointLevelIndex(HashGridPoint::level_type l, HashGridPoint::index_type i,
                                          level_t& lOut, HashGridPoint::index_type& iOut) {
  if (l == 0 && i == 0) {
    lOut = -1;
    iOut = 0;
  } else if (l == 0 && i == 1) {
    lOut = 0;
    iOut = 1;
  } else {
    lOut = l;
    iOut = i;
  }
}
  }  // namespace AnovaTypes
}  // namespace base
}  // namespace sgpp