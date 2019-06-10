#include "DimReduction.hpp"
#include <random>
#include "Tools.hpp"

namespace sgpp {
namespace base {


TransformationFunction::TransformationFunction(DataMatrix transformation) : VectorFunction(transformation.getNcols(), transformation.getNrows()), transformation(transformation) {
}

void TransformationFunction::eval(const base::DataVector& x, DataVector& out) {
  out = Tools::mult(transformation, x);
}

void TransformationFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::make_unique<TransformationFunction>(transformation);
}
}  // namespace base
}  // namespace sgpp
