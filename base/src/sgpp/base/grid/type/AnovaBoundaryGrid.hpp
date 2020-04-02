// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/AnovaBoundaryGridGenerator.hpp>
#include <sgpp/base/grid/AnovaTypes.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>

namespace sgpp {
namespace base {

/**
 * ANOVA grid with linear base functions with boundaries.
 */
class AnovaBoundaryGrid : public Grid {
 protected:
  explicit AnovaBoundaryGrid(std::istream& istr);

 public:
  /**
   * A false entry for a certain dimension indicates that the basis function is constant in that dimension.
   * 
   * A true entry for a certain dimension indicates that the basis function is active in that dimension.
   */
  typedef std::vector<bool> AnovaComponent;

  /**
   * A vector of ANOVA components.
   */
  typedef std::vector<AnovaComponent> AnovaComponentVector;

  
    /**
   * Gets the ANOVA component for a grid point.
   */
  static AnovaComponent getAnovaComponentOfPoint(const GridPoint& point)
{
    AnovaComponent currentComp(point.getDimension(), false);
    for (size_t d = 0; d < point.getDimension(); d++) {
      currentComp[d] = !(point.getLevel(d) == 0 && point.getIndex(d) == 0);
    }
    return currentComp;
  }

  /**
   * Level type for ANOVA grids (uses -1 level)
   */
  typedef int32_t level_t;

  /**
   * Converts a compatibility level to the true ANOVA grid level
   */
  static level_t fromNormalLevel(base::level_t l) {
    if (l == 0) {
      return -1;
    } else {
      return static_cast<level_t>(l - 1);
    }
  }

    /**
   * Converts an ANOVA grid level to the compatibility level
   */
  static base::level_t toNormalLevel(level_t l) {
    if (l == -1) {
      return 0;
    } else {
      return static_cast<base::level_t>(l + 1);
    }
  }

  /**
   * Converts the level and index ofa true ANOVA grid point to the level and index of a compatibility grid point
   */
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

  /**
   * Converts the level and index of a compatibility grid point to the level and index of a true
   * ANOVA grid point
   */
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
   * Constructor Anova Boundary Grid
   *
   * @param dim           the dimension of the grid
   */
  AnovaBoundaryGrid(size_t dim);

  /**
   * Destructor
   */
  ~AnovaBoundaryGrid() override = default;

  GridGenerator& getGenerator() override;

  void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION) override;

 protected:
  /// grid generator
  AnovaBoundaryGridGenerator generator;
};

}  // namespace base
}  // namespace sgpp
