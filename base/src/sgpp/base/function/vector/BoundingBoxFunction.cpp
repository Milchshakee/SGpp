#include <sgpp/base/function/vector/BoundingBoxFunction.hpp>

sgpp::base::BoundingBoxFunction::BoundingBoxFunction(Type type, const BoundingBox& bb)
    : VectorFunction(bb.getDimension(), bb.getDimension()), type(type), bb(bb) {}

void sgpp::base::BoundingBoxFunction::eval(const DataVector& x, DataVector& value) {
  if (type == Type::TO_UNIT_BB) {
    for (size_t d = 0; d < bb.getDimension(); d++) {
      value[d] = (x[d] - bb.getBoundary(d).leftBoundary) /
                 (bb.getBoundary(d).leftBoundary + bb.getBoundary(d).rightBoundary);
    }
  }

  else {
    for (size_t d = 0; d < bb.getDimension(); d++) {
      value[d] = (x[d] * (bb.getBoundary(d).leftBoundary +
                  bb.getBoundary(d).rightBoundary)) + bb.getBoundary(d).leftBoundary;
    }
  }
}

void sgpp::base::BoundingBoxFunction::clone(std::unique_ptr<VectorFunction>& clone) const {
  clone = std::make_unique<BoundingBoxFunction>(type, bb);
}
