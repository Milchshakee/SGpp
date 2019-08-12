#include <sgpp/base/grid/type/AnovaLinearBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/AnovaLinearBoundaryBasis.hpp>

sgpp::base::AnovaLinearBoundaryGrid::AnovaLinearBoundaryGrid(std::istream& istr) : AnovaBoundaryGrid(istr) {
}

sgpp::base::AnovaLinearBoundaryGrid::AnovaLinearBoundaryGrid(size_t dim)
    : AnovaBoundaryGrid(dim) {}

sgpp::base::GridType sgpp::base::AnovaLinearBoundaryGrid::getType() {
  return sgpp::base::GridType::AnovaLinearBoundary;
}

sgpp::base::SBasis& sgpp::base::AnovaLinearBoundaryGrid::getBasis() {
  static SAnovaLinearBoundaryBasis basis;
  return basis;
}

sgpp::base::Grid* sgpp::base::AnovaLinearBoundaryGrid::unserialize(std::istream& istr) {
  return new AnovaLinearBoundaryGrid(istr);
}
