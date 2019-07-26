#include "DimReduction.hpp"
#include <random>
#include "EigenHelper.hpp"

namespace sgpp {
namespace base {


TransformationFunction::TransformationFunction(DataMatrix transformation) : VectorFunction(transformation.getNcols(), transformation.getNrows()), transformation(transformation) {
}

void TransformationFunction::eval(const base::DataVector& x, DataVector& out) {
  out = EigenHelper::mult(transformation, x);
}

void TransformationFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::make_unique<TransformationFunction>(transformation);
}


size_t TransformationFunction::getOldDimensions() { return transformation.getNcols(); }

size_t TransformationFunction::getNewDimensions() { return transformation.getNrows(); }


EvalFunction::EvalFunction(const SGridSample& sample)
    : ScalarFunction(sample.getDimensions()), sample(&sample) {}

double EvalFunction::eval(const base::DataVector& x) { return sample->eval(x); }

void EvalFunction::clone(std::unique_ptr<ScalarFunction>& clone) const {
  clone = std::make_unique<EvalFunction>(*sample);
}
}  // namespace base
}  // namespace sgpp
