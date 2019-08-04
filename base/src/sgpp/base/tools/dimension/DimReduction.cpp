#include "DimReduction.hpp"
#include <random>
#include "EigenHelper.hpp"

namespace sgpp {
namespace base {

InputProjection::InputProjection(const DataMatrix& basis, size_t n, const DataVector& mean) : basis(basis), oldDimensions(basis.getNrows()), newDimensions(n), mean(mean) {
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


void InputProjection::inverse(DataVector& in, DataVector& out) {
    out = start;
  for (size_t d = 0; d < newDimensions; ++d) {
    double scale = posRange[d] - negRange[d];
    DataVector add(oldDimensions);
    basis.getColumn(d, add);
    add.mult(scale * in[d]);
    out.add(add);
  }
}

InputProjection::ProjectionFunction::ProjectionFunction(InputProjection& p) : VectorFunction(p.oldDimensions, p.newDimensions), p(p) {
}

void InputProjection::ProjectionFunction::eval(const DataVector& in, DataVector& out) {
}

void InputProjection::ProjectionFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
}
}  // namespace base
}  // namespace sgpp
