#include <random>
#include "MatrixFunction.hpp"
#include <sgpp/base/tools/dimension/EigenHelper.hpp>

namespace sgpp {
namespace base {

MatrixFunction::MatrixFunction(DataMatrix transformation)
    : VectorFunction(transformation.getNcols(), transformation.getNrows()),
      transformation(transformation) {}

void MatrixFunction::eval(const base::DataVector& x, DataVector& out) {
  out = EigenHelper::mult(transformation, x);
}

void MatrixFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::make_unique<MatrixFunction>(transformation);
}

size_t MatrixFunction::getOldDimensions() { return transformation.getNcols(); }

size_t MatrixFunction::getNewDimensions() { return transformation.getNrows(); }
}  // namespace base
}  // namespace sgpp