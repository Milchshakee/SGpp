#include "AsReducer.hpp"
#include "Tools.hpp"

sgpp::base::AsReducedFunction::AsReducedFunction(
    std::unique_ptr<sgpp::optimization::ScalarFunction>&& function, DataMatrix transformation) :
    ScalarFunction(transformation.getNcols()),
    function(std::move(function)), transformation(transformation) {}

double sgpp::base::AsReducedFunction::eval(const base::DataVector& x) {
  DataVector y = Tools::mult(transformation, x);
  return function->eval(y);
}

void sgpp::base::AsReducedFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  std::unique_ptr<ScalarFunction> ptr;
  function->clone(ptr);
  clone = std::make_unique<AsReducedFunction>(std::move(ptr), transformation);
}

sgpp::base::AsResult::AsResult(const DataMatrix& m, size_t n) {
  transformation = m;
  transformation.resizeRowsCols(m.getNrows(), n);
}

std::unique_ptr<sgpp::base::AsReducedFunction> sgpp::base::AsResult::apply(
  sgpp::optimization::ScalarFunction& input) {
  std::unique_ptr<optimization::ScalarFunction> ptr;
  input.clone(ptr);
  return std::make_unique<AsReducedFunction>(std::move(ptr), transformation);
}