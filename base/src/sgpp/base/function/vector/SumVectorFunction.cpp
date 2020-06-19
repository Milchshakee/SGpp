#include <sgpp/base/function/vector/SumVectorFunction.hpp>

sgpp::base::SumVectorFunction::SumVectorFunction(
    std::vector<std::shared_ptr<VectorFunction>>& functions,
                                     std::vector<bool> signs)
    : VectorFunction(functions.size() > 0 ? functions[0]->getNumberOfParameters() : 0,
                     functions.size() > 0 ? functions[0]->getNumberOfComponents() : 0),
      functions(functions),
      signs(signs) {}

void sgpp::base::SumVectorFunction::eval(const DataVector& x, DataVector& out) {
  for (size_t i = 0; i < functions.size(); i++)
  {
    double sign = signs.size() > 0 ? (signs[i] ? 1 : -1) : 1;
    DataVector r(m);
    functions[i]->eval(x, r);
    r.mult(sign);
    out.add(r);
  } 
}

void sgpp::base::SumVectorFunction::clone(std::unique_ptr<VectorFunction>& clone) const {}
