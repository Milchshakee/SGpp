#include <sgpp/base/function/scalar/EvalFunction.hpp>

namespace sgpp {
namespace base {


EvalFunction::EvalFunction() : ScalarFunction(0), sample(nullptr) {}

EvalFunction::EvalFunction(const SGridSample& sample)
    : ScalarFunction(sample.getDimensions()), sample(&sample)
{
  if (!sample.isHierarchised())
  {
    throw std::invalid_argument("Sample must be hierachised");
  } }

double EvalFunction::eval(const base::DataVector& x) { return sample->eval(x); }

void EvalFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::make_unique<EvalFunction>(*sample);
}
}  // namespace base
}  // namespace sgpp