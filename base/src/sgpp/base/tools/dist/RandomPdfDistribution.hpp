// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef RANDOMPDFDISTRIBUTION_HPP
#define RANDOMPDFDISTRIBUTION_HPP

#include <random>
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/tools/VectorDistribution.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {

class RandomPdfDistribution : public RandomDistribution {
 public:
  RandomPdfDistribution(size_t size, uint64_t seed, size_t dimensions, GridSample<double> pdf);

  void generate();

 private:
  struct Block
  {
    
  };

  size_t dimensions;
  std::vector<double> probs;
  std::vector<Block> blocks;
};

}  // namespace base
}  // namespace sgpp

#endif