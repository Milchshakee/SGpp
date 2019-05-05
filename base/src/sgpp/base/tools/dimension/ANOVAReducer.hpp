// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ANOVAREDUCER_HPP
#define ANOVAREDUCER_HPP

#include <sgpp/base/tools/dimension/DimReduction.hpp>

class ANOVAReducer : public sgpp::base::FunctionReducer {
 public:
  ANOVAReducer(size_t anovaOrder);
  std::unique_ptr<sgpp::optimization::ScalarFunction> reduceFunction(
      sgpp::optimization::ScalarFunction& input) override;

 private:
  size_t level;
  size_t anovaOrder;
};

#endif
