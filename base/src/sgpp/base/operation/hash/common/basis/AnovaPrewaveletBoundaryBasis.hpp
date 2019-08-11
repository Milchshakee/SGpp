// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <cmath>

namespace sgpp {
namespace base {

/**
 * Prewavelet basis on ANOVA boundary grids.
 */
template <class LT, class IT>
class AnovaPrewaveletBoundaryBasis : public Basis<LT, IT> {
 private:
  static const double normal_stamp[];
  static const double border_stamp[];

  double inline evalNormalHat(LT level, IT index, double p) {
    return 1.0 - fabs(static_cast<double>(1 << level) * p - static_cast<double>(index));
  }

 public:
  /**
   * Constructor
   */
  explicit AnovaPrewaveletBoundaryBasis() {}

  /**
   * Destructor.
   */
  ~AnovaPrewaveletBoundaryBasis() override {}

  /**
   * Evaluate a basis function.
   * Do not call this function for level -1, since it is constant at level -1.
   *
   * @param level     level of basis function
   * @param index     index of basis function
   * @param p     evaluation point
   * @return      value of Prewavelet basis function
   */
  inline double eval(LT level, IT index, double p) override {
    if (level == 0) {
      return 2 * (p - 0.5);
    }
    
    if (p == 0.0 || p == 1.0) {
      return 0.0;
    }

    if (1 == level) {
      return evalNormalHat(level, index, p);
    } else if (1 == index) {  // left border
      // Index of the affected hatbasis. The affected bases are ab and ab + 1
      int ab = static_cast<int>(floor(p * static_cast<double>(1 << level)));

      if (ab == 0) {
        return 0.9 * evalNormalHat(level, 1, p);
      } else if (ab == 3) {
        return 0.1 * evalNormalHat(level, 3, p);
      } else {
        return border_stamp[ab - 1] * evalNormalHat(level, ab, p) +
               border_stamp[ab] * evalNormalHat(level, ab + 1, p);
      }
    } else if ((unsigned int)(1 << level) - 1 == index) {  // right border
      // Index of the affected hatbasis. The affected bases are ab and ab + 1
      int ab = static_cast<int>(floor(p * static_cast<double>(1 << level)));

      if (ab == (1 << level) - 1) {
        return 0.9 * evalNormalHat(level, (1 << level) - 1, p);
      } else if (ab == (1 << level) - 4) {
        return 0.1 * evalNormalHat(level, (1 << level) - 3, p);
      } else {
        return border_stamp[(1 << level) - 1 - ab] * evalNormalHat(level, ab, p) +
               border_stamp[(1 << level) - 2 - ab] * evalNormalHat(level, ab + 1, p);
      }
    } else {
      // Index of the affected hatbasis. The affected bases are ab and ab + 1
      unsigned int ab = static_cast<int>(floor(p * static_cast<double>(1 << level)));

      if (ab == index - 3) {
        return normal_stamp[0] * evalNormalHat(level, ab + 1, p);
      } else if (ab == index + 2) {
        return normal_stamp[4] * evalNormalHat(level, ab, p);
      } else {
        int stamp = ab - index + 2;
        return normal_stamp[stamp] * evalNormalHat(level, ab, p) +
               normal_stamp[stamp + 1] * evalNormalHat(level, ab + 1, p);
      }
    }
  }

  double getIntegral(LT level, IT index) override { return -1.0;
  }

  size_t getDegree() const override { return 1; }
};

template <class LT, class IT>
const double AnovaPrewaveletBoundaryBasis<LT, IT>::normal_stamp[] = {0.1, -0.6, 1.0, -0.6, 0.1};

template <class LT, class IT>
const double AnovaPrewaveletBoundaryBasis<LT, IT>::border_stamp[] = {0.9, -0.6, 0.1};


// default type-def (unsigned int for level and index)
// This works because we don't need the level -1 (Because the value is always 1)
typedef AnovaPrewaveletBoundaryBasis<unsigned int, unsigned int> SAnovaPrewaveletBoundaryBasis;

}  // namespace base
}  // namespace sgpp
