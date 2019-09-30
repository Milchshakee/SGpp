#include <sgpp/base/tools/EigenHelper.hpp>
#include <sgpp/base/tools/dimension/DimReduction.hpp>

namespace sgpp {
namespace base {

InputProjection::InputProjection(const DataMatrix& basis, size_t n, const DataVector& mean)
    : oldDimensions(basis.getNrows()), newDimensions(n), mean(mean), func(*this) {
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
  posRange = DataVector(dim, 0);
  negRange = DataVector(dim, 0);
  for (size_t d = 0; d < dim; ++d) {
    for (size_t j = 0; j < oldDimensions; ++j) {
      DataVector col(oldDimensions);
      oldToNewBasis.getRow(d, col);
      double v = col[j];
      double m = mean[j];

      posRange[d] += v >= 0 ? (1 - m) * v : -m * v;
      negRange[d] += v < 0 ? (1 - m) * v : -m * v;
    }
  }

  start = mean;
  end = mean;
  for (size_t d = 0; d < newDimensions; ++d) {
    DataVector add(oldDimensions);
    oldToNewBasis.getRow(d, add);
    add.mult(negRange[d]);
    start.add(add);
    oldToNewBasis.getRow(d, add);
    add.mult(posRange[d]);
    end.add(add);
  }
}

void InputProjection::inverse(const DataVector& in, DataVector& out) {
  out = start;
  for (size_t d = 0; d < newDimensions; ++d) {
    double scale = posRange[d] - negRange[d];
    DataVector add(oldDimensions);
    oldToNewBasis.getRow(d, add);
    add.mult(scale * in[d]);
    out.add(add);
  }
}

const DataVector& InputProjection::getStart() { return start; }

const DataVector& InputProjection::getEnd() { return end; }

InputProjection::ProjectionFunction::ProjectionFunction(InputProjection& p)
    : VectorFunction(p.oldDimensions, p.newDimensions), p(&p) {}

void InputProjection::ProjectionFunction::eval(const DataVector& in, DataVector& out) {
  out = in;
  out.sub(p->mean);
  out = EigenHelper::mult(p->oldToNewBasis, out);
  DataVector ranges(p->newDimensions);
  for (size_t d = 0; d < p->newDimensions; ++d) {
    double scale = p->posRange[d] - p->negRange[d];
    if (out[d] > p->posRange[d]) {
      out[d] = p->posRange[d];
    } else if (out[d] < p->negRange[d]) {
      out[d] = p->negRange[d];
    }

    out[d] = (out[d] - p->negRange[d]) / scale;
  }
}

void InputProjection::ProjectionFunction::clone(std::unique_ptr<VectorFunction>& clone) const {}

VectorFunction& InputProjection::getFunction() { return func; }

L2SquaredMcRule::L2SquaredMcRule(uint64_t seed, size_t samples) : seed(seed), samples(samples) {}

L2McRule::L2McRule(uint64_t seed, size_t samples) : seed(seed), samples(samples) {}

double L2McRule::calculateRelativeError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r) {
  double e = calculateAbsoluteError(f, t, r);
  double fe = calculateAbsoluteError(f);
  return fe != 0.0 ? e / fe : 0.0;
}

double L2McRule::calculateAbsoluteError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r) {
  OperationL2 o(seed, samples);
  double e = o.calculateMcL2Error(f, t, r);
  return e;
}

double L2McRule::calculateAbsoluteError(ScalarFunction& f) {
  OperationL2 o(seed, samples);
  double l2 = o.calculateMcL2Norm(f);
  return l2;
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
