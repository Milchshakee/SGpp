// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <cmath>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>

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
    return std::max(1.0 - fabs(static_cast<double>(1 << level) * p - static_cast<double>(index)),
                    0.0);
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

    if (1 == level) {
      if (p <= 0.5) {
        return -1 + (p * 4);
      } else {
        return 1 - ((2 * p - 1) * 2);
      }
    }

    if (1 == index) {
      double res = 0;
      for (size_t i = 0; i < 4; i++) {
        res += border_stamp[i] * evalNormalHat(level, i, p);
      }
      return res;
    }

    if ((unsigned int)(1 << level) - 1 == index) {  // right border
      double res = 0;
      for (int32_t i = 3; i >= 0; i--) {
        res += border_stamp[i] * evalNormalHat(level, (1 << level) - i, p);
      }
      return res;
    }

    if (p == 0.0 || p == 1.0) {
      return 0.0;
    }

    // Index of the affected hatbasis. The affected bases are ab and ab + 1
    unsigned int ab = static_cast<int>(floor(p * static_cast<double>(1 << level)));
    if (ab >= index - 3 && ab <= index + 2) {
      if (ab == index - 3) {
        return normal_stamp[0] * evalNormalHat(level, ab + 1, p);
      } else if (ab == index + 2) {
        return normal_stamp[4] * evalNormalHat(level, ab, p);
      } else {
        int stamp = ab - index + 2;
        return normal_stamp[stamp] * evalNormalHat(level, ab, p) +
               normal_stamp[stamp + 1] * evalNormalHat(level, ab + 1, p);
      }
    } else {
      return 0.0;
    }
  }

  double getIntegral(LT level, IT index) override { return -1.0; }

  size_t getDegree() const override { return 1; }
};

template <class LT, class IT>
const double AnovaPrewaveletBoundaryBasis<LT, IT>::normal_stamp[] = {0.1, -0.6, 1.0, -0.6, 0.1};

template <class LT, class IT>
const double AnovaPrewaveletBoundaryBasis<LT, IT>::border_stamp[] = {-1.2, 1.1, -0.6, 0.1};

// default type-def (unsigned int for level and index)
// This works because we don't need the level -1 (Because the value is always 1)
typedef AnovaPrewaveletBoundaryBasis<unsigned int, unsigned int> SAnovaPrewaveletBoundaryBasis;

}  // namespace base
}  // namespace sgpp
