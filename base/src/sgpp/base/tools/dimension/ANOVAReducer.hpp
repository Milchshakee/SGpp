// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ANOVAREDUCER_HPP
#define ANOVAREDUCER_HPP

#include <sgpp/base/tools/dimension/DimReduction.hpp>

struct ANOVAInformation
{
  
};

class ANOVAReducer : public sgpp::base::FunctionReducer<ANOVAInformation> {
 public:
  ANOVAReducer(size_t anovaOrder);

protected:
  void evaluateFunction(sgpp::optimization::ScalarFunction& input, ANOVAInformation& out) override;
  std::unique_ptr<sgpp::optimization::ScalarFunction> reduce(
      sgpp::optimization::ScalarFunction& input, size_t n, const ANOVAInformation& info) override;

 private:
  size_t level;
  size_t anovaOrder;
};

#endif
