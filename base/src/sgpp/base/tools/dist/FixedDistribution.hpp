// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/dist/VectorDistribution.hpp>

namespace sgpp {
namespace base {

class FixedDistribution : public VectorDistribution {
 public:
  FixedDistribution(const DataMatrix& data);
  FixedDistribution(size_t size, size_t dimensions, const std::vector<DataVector>& vectors);
};
}  // namespace base
}  // namespace sgpp
