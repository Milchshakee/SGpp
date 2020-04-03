// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <random>
#include <vector>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>

namespace sgpp {
namespace base {

class DistributionSample {
 public:
  DistributionSample(const std::vector<DataVector>& vectors);
  DistributionSample(const DataMatrix& data);
  DistributionSample(sgpp::base::Grid& grid);
  DistributionSample(size_t size, sgpp::base::DistributionsVector& dist);
  DistributionSample(size_t size, std::default_random_engine& gen, ScalarFunction& pdf,
                     size_t iterations, double stepSize);

  size_t getSize() const;
  size_t getDimensions() const;
  const std::vector<DataVector>& getVectors() const;
  DataMatrix getAsDataMatrix() const;

 private:
  size_t size;
  size_t dimensions;
  std::vector<DataVector> vectors;
};
}  // namespace base
}  // namespace sgpp

