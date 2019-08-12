#include <sgpp/base/grid/type/AnovaPrewaveletBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/AnovaPrewaveletBoundaryBasis.hpp>

sgpp::base::AnovaPrewaveletBoundaryGrid::AnovaPrewaveletBoundaryGrid(std::istream& istr)
    : AnovaBoundaryGrid(istr) {}

sgpp::base::AnovaPrewaveletBoundaryGrid::AnovaPrewaveletBoundaryGrid(size_t dim)
    : AnovaBoundaryGrid(dim) {}

sgpp::base::GridType sgpp::base::AnovaPrewaveletBoundaryGrid::getType() {
  return sgpp::base::GridType::AnovaLinearBoundary;
}

sgpp::base::SBasis& sgpp::base::AnovaPrewaveletBoundaryGrid::getBasis() {
  static SAnovaPrewaveletBoundaryBasis basis;
  return basis;
}

sgpp::base::Grid* sgpp::base::AnovaPrewaveletBoundaryGrid::unserialize(std::istream& istr) {
  return new AnovaPrewaveletBoundaryGrid(istr);
}