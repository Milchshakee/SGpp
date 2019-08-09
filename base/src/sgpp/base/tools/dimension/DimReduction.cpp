#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/base/tools/dimension/AsQuadReducer.hpp>

namespace sgpp {
namespace base {

InputProjection::InputProjection(const DataMatrix& basis, size_t n, const DataVector& mean) : basis(basis), oldDimensions(basis.getNrows()), newDimensions(n), mean(mean), func(*this) {
  cutBasis = basis;
  cutBasis.resizeRowsCols(oldDimensions, newDimensions);
  cutBasis.transpose();
  calculateRanges();
}

  void InputProjection::calculateRanges() {
  size_t dim = newDimensions;
  posRange = DataVector(dim, 2);
  negRange = DataVector(dim, -2);
  for (size_t d = 0; d < dim; ++d) {
    double m = mean[d];
    for (size_t i = 0; i < oldDimensions; ++i) {
      DataVector col(oldDimensions);
      basis.getColumn(i, col);
      double v = col[d];
      double posX = (1 - m) / v;
      double negX = (0 - m) / v;
      if (posX < 0) {
        std::swap(posX, negX);
      }
      if (posX < posRange[d]) {
        posRange[d] = posX;
      }
      if (negX > negRange[d]) {
        negRange[d] = negX;
      }
    }
  }

    start = mean;
  for (size_t d = 0; d < newDimensions; ++d) {
    DataVector add(oldDimensions);
    basis.getColumn(d, add);
    add.mult(negRange[d]);
    start.add(add);
  }
  }


void InputProjection::inverse(const DataVector& in, DataVector& out) {
    out = start;
  for (size_t d = 0; d < newDimensions; ++d) {
    double scale = posRange[d] - negRange[d];
    DataVector add(oldDimensions);
    basis.getColumn(d, add);
    add.mult(scale * in[d]);
    out.add(add);
  }
}

InputProjection::ProjectionFunction::ProjectionFunction(InputProjection& p) : VectorFunction(p.oldDimensions, p.newDimensions), p(&p) {
}

void InputProjection::ProjectionFunction::eval(const DataVector& in, DataVector& out) {
  out = in;
  out.sub(p->mean);
  out = EigenHelper::mult(p->cutBasis, out);
  DataVector ranges(p->newDimensions);
  for (size_t d = 0; d < p->newDimensions; ++d) {
    double scale = p->posRange[d] - p->negRange[d];
    if (out[d] > p->posRange[d]) {
      out[d] = p->posRange[d];
      }
    else if (out[d] < p->negRange[d]) {
      out[d] = p->negRange[d];
    }

    out[d] = (out[d] - p->negRange[d]) / scale;
  }
}

void InputProjection::ProjectionFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
}


VectorFunction& InputProjection::getFunction() { return func; }


VarianceMcL2Rule::VarianceMcL2Rule(uint64_t seed, size_t samples) : seed(seed), samples(samples) {
}

double VarianceMcL2Rule::calculateRelativeError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r) {
  OperationL2 o(seed, samples);
  double e = o.calculateMcL2Error(f, t, r);
  double l2 = o.calculateMcL2Norm(f);
  return std::pow(e / l2, 2);
}


double VarianceMcL2Rule::calculateAbsoluteError(ScalarFunction& f, VectorFunction& t,
  ScalarFunction& r) {
  OperationL2 o(seed, samples);
  double e = o.calculateMcL2Error(f, t, r);
  return std::pow(e, 2);
}
}  // namespace base
}  // namespace sgpp
