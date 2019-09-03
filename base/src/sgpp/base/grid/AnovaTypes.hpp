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
}  // namespace AnovaTypes
}  // namespace base
}  // namespace sgpp