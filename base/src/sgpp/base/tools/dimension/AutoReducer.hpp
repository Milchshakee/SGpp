// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef AUTOREDUCER_HPP
#define AUTOREDUCER_HPP

#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include <sgpp/base/tools/dimension/Reducer.hpp>

class AutoReducer : public Reducer {
 public:
  std::unique_ptr<sgpp::optimization::VectorFunction> reduce(
      sgpp::optimization::VectorFunction& input);
};

#endif
