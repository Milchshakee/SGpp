// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ITERATIVEPCAREDUCER_HPP
#define ITERATIVEPCAREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "AsMcReducer.hpp"

class IterativePCAReducer : public sgpp::base::DataReducer {
 public:
  std::unique_ptr<sgpp::optimization::ScalarFunction> reduceData(sgpp::base::DataMatrix& input) override;

 private:
  size_t iterations;
  std::mt19937_64 prng;
};

#endif
