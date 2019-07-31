// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ANOVABOUNDARYGRID_HPP
#define ANOVABOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include "sgpp/base/grid/generation/AnovaBoundaryGridGenerator.hpp"

namespace sgpp {
namespace base {

/**
 * grid with linear base functions with boundaries, pentagon cut
 */
class AnovaBoundaryGrid : public Grid {
 public:
  typedef std::vector<bool> AnovaComponent;
  typedef std::vector<AnovaComponent> AnovaComponentVector;

  

static AnovaComponent getAnovaComponentOfPoint(const GridPoint& point)
{
    AnovaComponent currentComp(point.getDimension(), false);
    for (size_t d = 0; d < point.getDimension(); d++) {
      currentComp[d] = point.getLevel(d) > 0;
    }
    return currentComp;
  }

  typedef int32_t level_t;

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

  static level_t fromNormalLevel(base::level_t l) {
    if (l == 0) {
      return -1;
    } else {
      return static_cast<level_t>(l - 1);
    }
  }

  static base::level_t toNormalLevel(level_t l) {
    if (l == -1) {
      return 0;
    } else {
      return static_cast<base::level_t>(l + 1);
    }
  }

  static void toNormalGridPointLevelIndex(level_t l, HashGridPoint::index_type i,
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

  static void fromNormalGridPointLevelIndex(HashGridPoint::level_type l,
                                            HashGridPoint::index_type i, level_t& lOut,
                                            HashGridPoint::index_type& iOut) {
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

  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param dim           the dimension of the grid
   */
  AnovaBoundaryGrid(size_t dim);

  /**
   * Destructor
   */
  ~AnovaBoundaryGrid() override = default;

  sgpp::base::GridType getType() override;

  SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static Grid* unserialize(std::istream& istr);

  void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION) override;

 protected:
  /// grid generator
  AnovaBoundaryGridGenerator generator;
};

}  // namespace base
}  // namespace sgpp

#endif
