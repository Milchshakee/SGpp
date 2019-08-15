#include <sgpp/base/function/scalar/EvalFunction.hpp>

namespace sgpp {
namespace base {


EvalFunction::EvalFunction() : ScalarFunction(0), sample(nullptr) {}

EvalFunction::EvalFunction(const SGridSample& sample)
    : EvalFunction(sample, *op_factory::createOperationEval(const_cast<Grid&>(sample.getGrid()))) {
}

EvalFunction::EvalFunction(const SGridSample& sample, OperationEval& evalOp)
    : ScalarFunction(sample.getDimensions()), sample(&sample), evalOp(&evalOp)
{
  if (!sample.isHierarchised())
  {
    throw std::invalid_argument("Sample must be hierachised");
  } }

double EvalFunction::eval(const base::DataVector& x)
{ return evalOp->eval(sample->getValuesDataVector(), x); }

void EvalFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::make_unique<EvalFunction>(*sample, *evalOp);
}
}  // namespace base
}  // namespace sgpp