#pragma once

#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

class HypercubeFitFunction : public sgpp::base::VectorFunction
{

public:
  HypercubeFitFunction(sgpp::base::DataMatrix& basis, size_t reducedDims);
  ~HypercubeFitFunction() override;
  void eval(const sgpp::base::DataVector& x, sgpp::base::DataVector& value) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;

private:
  sgpp::base::DataMatrix basis;
 size_t reducedDims;

  sgpp::base::DataVector scale(sgpp::base::DataVector& point);

  sgpp::base::DataMatrix adjustBasis(sgpp::base::DataMatrix& basis, size_t reducedDims);

  sgpp::base::DataVector fromSphereCoords(const sgpp::base::DataVector& v);
  sgpp::base::DataVector toSphereCoords(const sgpp::base::DataVector& v);
};
