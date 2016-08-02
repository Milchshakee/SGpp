// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_

#include <sgpp/globaldef.hpp>
#include <functional>

namespace sgpp {
namespace combigrid {

class SingleFunction {
 public:
  typedef std::function<double(double)> function_type;

 private:
  function_type func;

 public:
  /**
   * for function pointers
   */
  SingleFunction(double (*ptr)(double));

  /**
   * for lambdas or function objects
   */
  template <typename T>
  explicit SingleFunction(T f)
      : func(f) {}

  double operator()(double param);
  double call(double param);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_ */