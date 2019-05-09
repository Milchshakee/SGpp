#include <iostream>

#include "sgpp/base/datatypes/DataMatrix.hpp"
#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include "sgpp/base/tools/dimension/ActiveSubspaceReducer.hpp"
#include "sgpp/optimization/function/vector/WrapperVectorFunction.hpp"
#include "sgpp/optimization/function/scalar/WrapperScalarFunction.hpp"

double f(const sgpp::base::DataVector& v)
{ return v[0]; }

void f_gradient(const sgpp::base::DataVector& x, sgpp::base::DataVector& out) {
  out[0] = 1;
  out[1] = 0;
}

int main() {
  auto func = std::make_shared<sgpp::optimization::WrapperScalarFunction>(2, f);
  auto funcGradient = std::make_shared<sgpp::optimization::WrapperVectorFunction>(2, 2, f_gradient);
  auto gradient = std::make_shared < sgpp::base::ActiveSubspaceReducer::GivenGradient>(funcGradient);
  auto dist = std::make_shared<sgpp::base::UniformVectorDistribution>(0, 2);
  auto cutoff = std::make_shared < sgpp::base::FixedCutoff<sgpp::base::ActiveSubspaceInfo>>(1);
  auto reducer = sgpp::base::ActiveSubspaceReducer(100, gradient, dist, cutoff);

  auto i = sgpp::base::ActiveSubspaceInfo();
  auto reducedFunc = reducer.reduceFunction(*func.get(), i);
  std::cout << reducedFunc->getNumberOfParameters() << std::endl;

  std::cout << reducedFunc->eval(sgpp::base::DataVector(1, 0.4)) << std::endl;
}
