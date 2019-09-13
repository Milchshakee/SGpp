// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>

namespace sgpp {
namespace base {

/**
 * ANOVA grid with prewavelet basis functions with boundaries.
 */
class AnovaPrewaveletBoundaryGrid : public AnovaBoundaryGrid {
 protected:
  explicit AnovaPrewaveletBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor Anova Boundary Grid
   *
   * @param dim           the dimension of the grid
   */
  AnovaPrewaveletBoundaryGrid(size_t dim);

    /**
   * Constructor Anova Boundary Grid
   *
   * @param dim           the dimension of the grid
   */
  AnovaPrewaveletBoundaryGrid(size_t dim, std::vector<AnovaTypes::LevelIndexPair>& anchor);

  /**
   * Destructor
   */
  ~AnovaPrewaveletBoundaryGrid() override = default;

  sgpp::base::GridType getType() override;

  SBasis& getBasis() override;

  static Grid* unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp
