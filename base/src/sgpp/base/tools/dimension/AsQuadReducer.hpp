// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASQUADREDUCER_HPP
#define ASQUADREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "AsReducer.hpp"
#include "sgpp/base/operation/hash/OperationQuadrature.hpp"

namespace sgpp {
namespace base {

    class AsQuadFixedCutter : public Cutter<GridSample<DataVector>, AsInfo, AsResult>
    {
    public:
  AsQuadFixedCutter(size_t n);

      AsResult cut(const GridSample<DataVector>& input, const AsInfo& info) override;

    private:
      size_t n;
    };

  class AsQuadEigenValueCutter : public AsEigenValueCutter<GridSample<DataVector>>
  {
  public:
    AsQuadEigenValueCutter(double minEigenValue);
  };

class AsQuadReducer : public AsReducer<GridSample<DataVector>> {
 public:
  AsInfo evaluate(GridSample<DataVector>& input) override;

};

}  // namespace base
}  // namespace sgpp

#endif
