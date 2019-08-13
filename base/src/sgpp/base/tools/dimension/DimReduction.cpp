#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/base/tools/EigenHelper.hpp>

namespace sgpp {
namespace base {

InputProjection::InputProjection(const DataMatrix& basis, size_t n, const DataVector& mean) : oldDimensions(basis.getNrows()), newDimensions(n), mean(mean), func(*this) {
  oldToNewBasis = basis;
  oldToNewBasis.transpose();
  oldToNewBasis.resizeRowsCols(newDimensions, oldDimensions);
  newToOldBasis = basis;
  newToOldBasis.transpose();
  newToOldBasis.resizeRowsCols(oldDimensions, newDimensions);
  calculateRanges();
}

  void InputProjection::calculateRanges() {
  size_t dim = newDimensions;
  posRange = DataVector(dim, 2);
  negRange = DataVector(dim, -2);
  for (size_t d = 0; d < dim; ++d) {
    double m = mean[d];
    for (size_t i = 0; i < newDimensions; ++i) {
      DataVector col(oldDimensions);
      newToOldBasis.getColumn(i, col);
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
    newToOldBasis.getColumn(d, add);
    add.mult(negRange[d]);
    start.add(add);
  }
  }


void InputProjection::inverse(const DataVector& in, DataVector& out) {
    out = start;
  for (size_t d = 0; d < newDimensions; ++d) {
    double scale = posRange[d] - negRange[d];
    DataVector add(oldDimensions);
    newToOldBasis.getColumn(d, add);
    add.mult(scale * in[d]);
    out.add(add);
  }
}

InputProjection::ProjectionFunction::ProjectionFunction(InputProjection& p) : VectorFunction(p.oldDimensions, p.newDimensions), p(&p) {
}

void InputProjection::ProjectionFunction::eval(const DataVector& in, DataVector& out) {
  out = in;
  out.sub(p->mean);
  out = EigenHelper::mult(p->oldToNewBasis, out);
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


L2SquaredMcRule::L2SquaredMcRule(uint64_t seed, size_t samples) : seed(seed), samples(samples) {
}


double ErrorRule::calculateRelativeError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r) {
  double e = calculateAbsoluteError(f, t, r);
  double fe = calculateAbsoluteError(f);
  return fe != 0.0 ? e / fe : 0.0;
}

double L2SquaredMcRule::calculateAbsoluteError(ScalarFunction& f, VectorFunction& t,
  ScalarFunction& r) {
  OperationL2 o(seed, samples);
  double e = o.calculateMcL2Error(f, t, r);
  return std::pow(e, 2);
}

  double L2SquaredMcRule::calculateAbsoluteError(ScalarFunction& f) {
  OperationL2 o(seed, samples);
  double l2 = o.calculateMcL2Norm(f);
  return std::pow(l2, 2);
}
}  // namespace base
}  // namespace sgpp
