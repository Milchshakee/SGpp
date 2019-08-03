// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASQUADREDUCER_HPP
#define ASQUADREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "AsReducer.hpp"
#include "sgpp/base/operation/hash/OperationQuadrature.hpp"
#include <sgpp/base/function/scalar/EvalFunction.hpp>

namespace sgpp {
namespace base {

  struct AsQuadInput
  {
  SGridSample functionSample;
  GridSample<DataVector> gradientSample;
  };

  class AsQuadResult : public AsResult<SGridSample>
  {
  public:
    AsQuadResult(const SGridSample& input, const DataMatrix& m, size_t n);

    ScalarFunction& getReducedFunction() override;
    SGridSample& getReducedOutput() override;
  private:
    SGridSample reduced;
   EvalFunction evalFunc;
  };

    class AsQuadFixedCutter : public Cutter<AsQuadInput, AsInfo, AsQuadResult>
    {
    public:
  AsQuadFixedCutter(size_t n);

      AsQuadResult cut(const AsQuadInput& input, const AsInfo& info) override;

    private:
      size_t n;
    };

class AsQuadReducer : public AsReducer<AsQuadInput> {
 public:
  AsInfo evaluate(AsQuadInput& input) override;

};

}  // namespace base
}  // namespace sgpp

#endif
