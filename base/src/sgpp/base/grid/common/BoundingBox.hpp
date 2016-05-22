// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BOUNDINGBOX_HPP
#define BOUNDINGBOX_HPP

#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

/**
 * struct that defines the boundaries for one specific dimension
 */
struct BoundingBox1D {
  /// left boundary
  double leftBoundary;
  /// right boundary
  double rightBoundary;
  /// whether to use Dirichlet boundaries on the left boundary
  bool bDirichletLeft;
  /// whether to use Dirichlet boundaries on the right boundary
  bool bDirichletRight;

  /**
   * Default constructor initializing leftBoundary = 0, rightBoundary = 1, and
   * bDirichletLeft = bDirichletLeft = false.
   */
  BoundingBox1D() : BoundingBox1D(0.0, 1.0, false, false) {}

  /**
   * Constructor initializing bDirichletLeft = bDirichletLeft = false.
   *
   * @param leftBoundary  left boundary position
   * @param rightBoundary right boundary position
   */
  BoundingBox1D(double leftBoundary, double rightBoundary) :
    BoundingBox1D(leftBoundary, rightBoundary, false, false) {}

  /**
   * Constructor.
   *
   * @param leftBoundary    left boundary position
   * @param rightBoundary   right boundary position
   * @param bDirichletLeft  whether to use Dirichlet boundaries on the left boundary
   * @param bDirichletRight whether to use Dirichlet boundaries on the right boundary
   */
  BoundingBox1D(double leftBoundary, double rightBoundary,
                bool bDirichletLeft, bool bDirichletRight) :
                  leftBoundary(leftBoundary), rightBoundary(rightBoundary),
                  bDirichletLeft(bDirichletLeft), bDirichletRight(bDirichletRight) {}
};

/**
 * This class implements the boundaries of the sparse grid.
 * Internally the grid is set up on a trivial cube.
 *
 * This class gives the class gives the opportunity to stretch
 * this cube in every dimension separately.
 */
class BoundingBox {
 protected:
  /// the number of dimensions used with the grid
  size_t dimension;
  /// Array that contains all left boundaries for all dimensions
  std::vector<BoundingBox1D> boundingBox1Ds;

 public:
  /**
   * Constructor for BoundingBox
   *
   * initializes the Bounding with a N-d trivial cube
   *
   * @param dimension number of the dimensions used with the grid
   */
  explicit BoundingBox(size_t dimension);

  /**
   * Constructor for BoundingBox
   *
   * initializes the Bounding with specific values for all dimensions
   *
   * @param boundingBox1Ds array that contains all boundaries
   */
  explicit BoundingBox(const std::vector<BoundingBox1D>& boundingBox1Ds);

  /**
   * Sets left and right boundary for a specific dimension.
   *
   * @param d             the dimension in which the boundary should be changed
   * @param boundingBox1D reference to a BoundingBox1D object that contains the new boundaries
   */
  void setBoundary(size_t d, const BoundingBox1D& boundingBox1D);

  /**
   * Returns the left and right boundary for a specific dimension.
   *
   * @param d   the dimension in which the boundary should be read
   * @return    a BoundingBox1D object that contains the boundaries
   */
  BoundingBox1D getBoundary(size_t d) const;

  /**
   * Returns the number of dimensions of this bounding box.
   *
   * @return number of dimensions
   */
  size_t getDimensions() const;

  /**
   * Calculates the width of the interval in one dimension.
   *
   * @param d the dimension in which the width of the interval should be determined
   * @return  width of the interval
   */
  double getIntervalWidth(size_t d) const;

  /**
   * Returns the offset in positive x-direction of the interval in one dimension.
   *
   * @param d dimension in which the offset of the interval should be determined
   * @return  offset in positive x-direction of the interval
   */
  double getIntervalOffset(size_t d) const;

  /**
   * Determine if this bounding box describes the unit cube \f$[0, 1]^d\f$.
   *
   * @return true if this bounding box is the unit cube, otherwise false
   */
  bool isUnitCube() const;

  /**
   * Determines if the interval in the specified dimension has left Dirichlet boundary conditions.
   *
   * @param dimension the dimension for which the left boundary condition should be determined
   * @return true     if Dirichlet Boundary conditions, otherwise false
   */
  bool hasDirichletBoundaryLeft(size_t d) const;

  /**
   * Determines if the interval in the specified dimension has right Dirichlet boundary conditions.
   *
   * @param dimension the dimension for which the right boundary condition should be determined
   * @return true     if Dirichlet Boundary conditions, otherwise false
   */
  bool hasDirichletBoundaryRight(size_t d) const;

  /**
   * Converts the BoundingBox to a string.
   *
   * @param text string to which the data is written
   */
  void toString(std::string& text) const;

  /**
   * Converts the BoundingBox to a string.
   *
   * @return string to which the data is written
   */
  std::string toString() const;
};

}  // namespace base
}  // namespace sgpp

#endif /* BOUNDINGBOX_HPP */
