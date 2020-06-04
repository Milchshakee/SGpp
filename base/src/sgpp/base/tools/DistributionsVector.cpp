// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>
#include <vector>
#include <sgpp/base/tools/DistributionTruncNormal.hpp>
#include <sgpp/base/function/vector/BoundingBoxFunction.hpp>
#include <sgpp/base/tools/DistributionTruncLogNormal.hpp>

namespace sgpp {
namespace base {

DistributionsVector::DistributionsVector() : distributions(0) {}

DistributionsVector::DistributionsVector(size_t dim) : distributions(dim) {}

DistributionsVector::DistributionsVector(size_t dim, std::shared_ptr<sgpp::base::Distribution> pdf)
    : distributions(dim, pdf) {}

DistributionsVector::DistributionsVector(const DistributionsVector& other) {
  distributions.clear();
  for (auto& pdf : other.distributions) {
    distributions.push_back(pdf);
  }
  bb = other.bb;
}

DistributionsVector::DistributionsVector(std::vector<DistributionType> types,
                                         std::vector<DataVector> characteristics,
                                         std::shared_ptr<BoundingBox> bb)
    : distributions(types.size()), bb(bb) {
  for (size_t d = 0; d < bb->getDimension(); d++) {
    double left = bb->getBoundary(d).leftBoundary;
    double right = bb->getBoundary(d).rightBoundary;
    std::shared_ptr<Distribution> dist;
    if (types[d] == DistributionType::Uniform) {
      dist = std::make_shared<DistributionUniform>(left, right);
    }
    else if (types[d] == DistributionType::TruncNormal || types[d] == DistributionType::TruncLognormal) {
      double mean = characteristics[d][0];
      double stddev = characteristics[d][1];
      if (types[d] == DistributionType::TruncNormal) {
        dist = std::make_shared<DistributionTruncNormal>(mean, stddev, left, right);
      } else {
        dist = std::make_shared<DistributionTruncLogNormal>(mean, stddev, left, right);
        }
    } else {
      throw std::invalid_argument("Invalid distribution");
      }
    distributions[d] = dist;
  }
}

DistributionsVector::~DistributionsVector() {}

std::vector<std::shared_ptr<sgpp::base::Distribution>> DistributionsVector::getDistributions() {
  return distributions;
}

sgpp::base::DataVector DistributionsVector::sample() {
  sgpp::base::DataVector sampleVector(distributions.size(), 0.0);
  for (size_t d = 0; d < distributions.size(); d++) {
    double v = distributions[d]->sample();
    if (bb != nullptr) {
    v = (v - bb->getBoundary(d).leftBoundary) /
               (bb->getBoundary(d).rightBoundary - bb->getBoundary(d).leftBoundary);
      }
    sampleVector[d] = v;
  }
  return sampleVector;
}

void DistributionsVector::push_back(std::shared_ptr<sgpp::base::Distribution> pdf) {
  distributions.push_back(pdf);
}

std::shared_ptr<sgpp::base::Distribution> DistributionsVector::get(size_t i) {
  return distributions[i];
}

size_t DistributionsVector::getSize() { return distributions.size(); }

void DistributionsVector::clear() { distributions.clear(); }

}  // namespace base
} /* namespace sgpp */
