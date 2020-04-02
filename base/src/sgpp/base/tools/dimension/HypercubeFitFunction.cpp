#include "HypercubeFitFunction.hpp"
#include <sgpp/base/tools/EigenHelper.hpp>

using namespace sgpp::base;

HypercubeFitFunction::HypercubeFitFunction(sgpp::base::DataMatrix& basis, size_t reducedDims)
    : VectorFunction(basis.getNcols(), reducedDims),
      basis(basis),
      reducedDims(reducedDims) {}

HypercubeFitFunction::~HypercubeFitFunction() {}

void HypercubeFitFunction::eval(const sgpp::base::DataVector& x, sgpp::base::DataVector& value)
{
  DataVector sphere = toSphereCoords(x);
  DataVector rotated = EigenHelper::mult(basis, sphere);
  value = fromSphereCoords(rotated);
}

void HypercubeFitFunction::clone(std::unique_ptr<VectorFunction>& clone) const {}

sgpp::base::DataMatrix HypercubeFitFunction::adjustBasis(sgpp::base::DataMatrix& basis,
                                                         size_t reducedDims) {
  DataMatrix transposed = basis;
  transposed.transpose();
  DataMatrix unitVectors = EigenHelper::mult(basis, transposed);

  size_t dims = basis.getNcols();
  DataMatrix cut = basis;
  cut.resizeRowsCols(cut.getNrows(), reducedDims);
  cut.transpose();
  //for (size_t i = reducedDims; i < dims; i++) {
  //  cut.setColumn(i, DataVector(dims, 0.0));
  //}

  DataMatrix newBasis = basis;
  for (size_t i = 0; i < reducedDims; i++) {
    DataVector unit(dims);
    basis.getColumn(i, unit);
    unit = EigenHelper::mult(transposed, unit);

    for (size_t j = reducedDims; j < dims; j++) {
      DataVector col(dims);
      basis.getColumn(j, col);
      DataVector toAdd(dims);
      double dot = unit.dotProduct(col);
      col.mult(dot);
      unit.add(col);
    }

    double a = 0;
    //newBasis.setColumn(i, projected);
  }
  return newBasis;
}


sgpp::base::DataVector HypercubeFitFunction::scale(sgpp::base::DataVector& point) {
  size_t dims = basis.getNcols();
  sgpp::base::DataVector base = point;
  for (size_t i = reducedDims; i < dims; i++) {
    base[i] = 0.5;
  }
  DataVector relative = point;
  relative.sub(base);
  relative.normalize();
  relative.mult(0.5);
  }


sgpp::base::DataVector HypercubeFitFunction::toSphereCoords(const sgpp::base::DataVector& v) {
    size_t dims = basis.getNcols();
    sgpp::base::DataVector base(dims, 0.5);
    DataVector relative = v;
    relative.sub(base);
    double sum = relative.l2Norm();
    relative.mult(1.0 / sum);
    return relative;
  }


sgpp::base::DataVector HypercubeFitFunction::fromSphereCoords(const sgpp::base::DataVector& v) {
    double sum = v.l2Norm();
  DataVector copy = v;
    copy.mult(sum);
    return copy;
  }
