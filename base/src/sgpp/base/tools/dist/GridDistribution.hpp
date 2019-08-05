// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <random>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <vector>
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "VectorDistribution.hpp"

namespace sgpp {
namespace base {

class GridDistribution : public VectorDistribution {
 public:
  GridDistribution(sgpp::base::Grid& grid);
};
}  // namespace base
}  // namespace sgpp
