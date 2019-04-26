// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REDUCER_HPP
#define REDUCER_HPP

#include "sgpp/optimization/function/vector/VectorFunction.hpp"

class Reducer
{
public:
  virtual ~Reducer();
  virtual std::unique_ptr<sgpp::optimization::VectorFunction> reduce(sgpp::optimization::VectorFunction& input);
};

#endif
