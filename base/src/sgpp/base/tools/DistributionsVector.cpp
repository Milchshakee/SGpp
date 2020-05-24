// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>
#include <vector>

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
}

DistributionsVector::DistributionsVector(std::vector<DistributionType> types, BoundingBox& bb,
                                         bool toUnitBB)
    : distributions(bb.getDimension()) {
  for (size_t d = 0; d < bb.getDimension(); d++) {
    std::shared_ptr<Distribution> dist;
    if (types[d] == DistributionType::Uniform) {
      dist = std::make_shared<DistributionUniform>(toUnitBB ? 0 : bb.getBoundary(d).leftBoundary,
                                                    toUnitBB ? 1 : bb.getBoundary(d).rightBoundary);
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
    sampleVector[d] = distributions[d]->sample();
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
