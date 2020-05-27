// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/Distribution.hpp>

#include <iostream>
#include <random>

namespace sgpp {
namespace base {

/**
 * Truncated LogNormal distribution.
 * Only accepts samples within [lower,upper]
 */
class DistributionTruncLogNormal : public Distribution {
 public:
  /**
   * Constructor
   */
  DistributionTruncLogNormal(double mean, double stddev, double lower, double upper)
      : Distribution(),
        mean(mean),
        stddev(stddev),
        lower(lower),
        upper(upper),
        dist(mean, stddev) {}

  /**
   * Destructor
   */
  virtual ~DistributionTruncLogNormal() {}

  /**
   *
   */
  double sample() {
    while (true) {
      double number = dist(gen);
      if (number >= lower && number <= upper) return number;
    }
  }

  /**
   *
   */
  double eval(double x) {
    if (x < lower) {
      std::cout << "DistributionTruncLogNormal: argument not in specified interval!\n";
      return lower;
    } else if (x > upper) {
      std::cout << "DistributionTruncLogNormal: argument not in specified interval!\n";
      return upper;
    } else {
      if (x >= 0)
        // y =   1.0./(sigma2*x*sqrt(2*pi)) .* exp(-((log(x)-mu).^2)./(2*sigma2^2));
        return 1.0 / (stddev * x * sqrt(2 * M_PI)) *
               exp(-std::pow(log(x) - mean, 2) / (2 * stddev * stddev));
      else
        return 0.0;
    }
  }

  /**
   *
   */
  sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds(2);
    bounds[0] = lower;
    bounds[1] = upper;
    return bounds;
  }

  sgpp::base::DistributionType getType() { return sgpp::base::DistributionType::TruncLognormal; }

  sgpp::base::DataVector getCharacteristics() {
    sgpp::base::DataVector characteristics(2);
    characteristics[0] = mean;
    characteristics[1] = stddev;
    return characteristics;
  }

 private:
  // mean, often called mu
  double mean;
  // standard deviation, often called sigma
  // Note: This is sigma, not sigma^2!
  double stddev;
  // lower bound
  double lower;
  // upper bound
  double upper;
  std::lognormal_distribution<double> dist;
};  // namespace base

}  // namespace base
}  // namespace sgpp
