// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/AnovaBoundaryGridGenerator.hpp>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/AnovaLinearBoundaryBasis.hpp>

namespace sgpp {
namespace base {

/**
 * ANOVA grid with linear base functions with boundaries.
 */
class AnovaLinearBoundaryGrid : public AnovaBoundaryGrid {
 protected:
  explicit AnovaLinearBoundaryGrid(std::istream& istr);

 public:

  /**
   * Constructor Anova Boundary Grid
   *
   * @param dim           the dimension of the grid
   */
  AnovaLinearBoundaryGrid(size_t dim);

  /**
   * Destructor
   */
  ~AnovaLinearBoundaryGrid() override = default;

  sgpp::base::GridType getType() override;

  SBasis& getBasis() override;

  static Grid* unserialize(std::istream& istr);

};

}  // namespace base
}  // namespace sgpp
