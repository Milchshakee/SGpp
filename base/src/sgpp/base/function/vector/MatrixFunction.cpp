#include <sgpp/base/tools/EigenHelper.hpp>
#include <sgpp/base/function/vector/MatrixFunction.hpp>

namespace sgpp {
namespace base {

MatrixFunction::MatrixFunction() : VectorFunction(0, 0) {}

MatrixFunction::MatrixFunction(DataMatrix transformation)
    : VectorFunction(transformation.getNcols(), transformation.getNrows()),
      transformation(transformation) {}

void MatrixFunction::eval(const base::DataVector& x, DataVector& out) {
  out = EigenHelper::mult(transformation, x);
}

void MatrixFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::make_unique<MatrixFunction>(transformation);
}

}  // namespace base
}  // namespace sgpp