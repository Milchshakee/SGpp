// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/Distribution.hpp>

#include <vector>
#include <sgpp/base/grid/common/BoundingBox.hpp>

namespace sgpp {
namespace base {

class DistributionsVector {
 public:
  DistributionsVector();
  explicit DistributionsVector(size_t dim);
  DistributionsVector(size_t dim, std::shared_ptr<sgpp::base::Distribution> pdf);
  DistributionsVector(const DistributionsVector& other);
  DistributionsVector(std::vector<DistributionType> types, BoundingBox& bb, bool toUnitBB);

  virtual ~DistributionsVector();

  std::vector<std::shared_ptr<sgpp::base::Distribution>> getDistributions();

  sgpp::base::DataVector sample();

  void push_back(std::shared_ptr<sgpp::base::Distribution> pdf);
  std::shared_ptr<sgpp::base::Distribution> get(size_t i);
  size_t getSize();
  void clear();

 private:
  std::vector<std::shared_ptr<sgpp::base::Distribution>> distributions;
};

}  // namespace base
} /* namespace sgpp */
