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
 protected:
  explicit AnovaBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param dim           the dimension of the grid
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  explicit AnovaBoundaryGrid(size_t dim);

  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param BB the BoundingBox of the grid
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  explicit AnovaBoundaryGrid(BoundingBox& BB);

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
