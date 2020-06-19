#include <sgpp/base/function/scalar/SumScalarFunction.hpp>

sgpp::base::SumScalarFunction::SumScalarFunction(
    std::vector<std::shared_ptr<ScalarFunction>>& functions,
                                     std::vector<bool> signs)
    : ScalarFunction(functions.size() > 0 ? functions[0]->getNumberOfParameters() : 0), functions(functions), signs(signs) {}

double sgpp::base::SumScalarFunction::eval(const DataVector& x) {
  double sum = 0;
  for (size_t i = 0; i < functions.size(); i++)
  {
    double sign = signs.size() > 0 ? (signs[i] ? 1 : -1) : 1;
    sum += sign * functions[i]->eval(x);
  } 
return sum;}

void sgpp::base::SumScalarFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {}
