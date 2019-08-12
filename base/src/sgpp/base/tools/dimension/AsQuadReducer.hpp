// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASQUADREDUCER_HPP
#define ASQUADREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/EvalFunction.hpp>
#include <sgpp/base/tools/dimension/AsReducer.hpp>

namespace sgpp {
namespace base {

  struct AsQuadInput
  {
  SGridSample functionSample;
  GridSample<DataVector> gradientSample;
  };

  class AsQuadResult : public Result<SGridSample>
  {
  public:
    AsQuadResult(const SGridSample& input, const DataMatrix& m, size_t n);


    ScalarFunction& getOriginalFunction() override;
    VectorFunction& getTransformationFunction() override;
    ScalarFunction& getReducedFunction() override;
    SGridSample& getReducedOutput() override;
  private:
    InputProjection projection;
    SGridSample reduced;
   EvalFunction evalFunc;
    EvalFunction originalFunc;
  };

    class AsQuadFixedCutter : public FixedCutter<AsQuadInput, AsInfo, AsQuadResult>
    {
    public:
  AsQuadFixedCutter(size_t n);

      AsQuadResult cut(const AsQuadInput& input, const AsInfo& info) override;
    };

class AsQuadReducer : public Reducer<AsQuadInput, AsInfo, AsQuadResult> {
 public:
  AsInfo evaluate(AsQuadInput& input) override;

};

}  // namespace base
}  // namespace sgpp

#endif
