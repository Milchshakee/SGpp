#include <sgpp/base/function/scalar/ChainScalarFunction.hpp>

sgpp::base::ChainScalarFunction::ChainScalarFunction(
    const std::vector<std::shared_ptr<VectorFunction>>& chain, const std::shared_ptr<ScalarFunction>& func)
    : ScalarFunction(chain[0]->getNumberOfParameters()), chain(chain), func(func) {}

sgpp::base::ChainScalarFunction::~ChainScalarFunction() {}

double sgpp::base::ChainScalarFunction::eval(const DataVector& x) {
  DataVector param = x;
  DataVector value;
  for (std::shared_ptr<VectorFunction>& v : chain) {
    value = DataVector(v->getNumberOfParameters());
    v->eval(param, value);
    param = value;
  }
  return func->eval(param);
}

void sgpp::base::ChainScalarFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::make_unique<ChainScalarFunction>(chain, func);
}
