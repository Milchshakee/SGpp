// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cmath>

namespace sgpp {
namespace base {

/**
 * Linear basis on ANOVA boundary grids.
 */
template <class LT, class IT>
class AnovaLinearBoundaryBasis : public Basis<LT, IT> {
 public:
  /**
   * Constructor
   */
  explicit AnovaLinearBoundaryBasis() {}

  /**
   * Destructor.
   */
  ~AnovaLinearBoundaryBasis() override {}

    double evalDx(LT level, IT index, double x) override {
    throw std::runtime_error("Not implemented");
  }

  /**
   * Evaluate a basis function.
   * Do not call this function for level -1, since it is constant at level -1.
   *
   * @param l     level of basis function
   * @param i     index of basis function
   * @param x     evaluation point
   * @return      value of linear basis function
   */
  inline double eval(LT l, IT i, double x) override {
    if (l == 0) {
      return x;
    } else {
      return std::max(1.0 - std::abs(static_cast<double>(static_cast<IT>(1) << (l)) * x -
                                     static_cast<double>(i)),
                      0.0);
    }
  }
  /**
   * Do not call this function for level -1
   */
  double getIntegral(LT level, IT index) override {
    if (level == 0) {
      return 0.5;
    } else {
      return 1. / static_cast<double>(static_cast<IT>(1) << (level));
    }
  }

  size_t getDegree() const override { return 1; }
};

// default type-def (unsigned int for level and index)
// This works because we don't need the level -1 (Because the value is always 1)
typedef AnovaLinearBoundaryBasis<unsigned int, unsigned int> SAnovaLinearBoundaryBasis;

}  // namespace base
}  // namespace sgpp
