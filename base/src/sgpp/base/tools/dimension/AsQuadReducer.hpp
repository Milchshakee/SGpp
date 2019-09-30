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
  ScalarFunction& originalFunction;
  SGridSample functionSample;
  GridSample<DataVector> gradientSample;
  };

  class AsQuadResult : public Result<SGridSample>
  {
  public:
    AsQuadResult(const AsQuadInput& input, const DataMatrix& m, size_t n, const DataVector& mean);

    
    ScalarFunction& getOriginalFunction() override;
    VectorFunction& getTransformationFunction() override;
    ScalarFunction& getReducedFunctionSurrogate() override;
    SGridSample& getReducedOutput() override;
    InputProjection& getProjection();

  private:
    ScalarFunction* originalFunction;
    InputProjection projection;
    SGridSample reduced;
   EvalFunction evalFunc;
  };

    class AsQuadFixedCutter : public FixedCutter<AsQuadInput, AsInfo, AsQuadResult>
    {
    public:
      AsQuadFixedCutter(size_t n, const DataVector& mean);

      AsQuadResult cut(const AsQuadInput& input, const AsInfo& info) override;

    private:
      DataVector mean;
    };

class AsQuadReducer : public Reducer<AsQuadInput, AsInfo, AsQuadResult> {
 public:
  static GridSample<DataVector> fromGradientFunction(std::shared_ptr<Grid>&,
                                                     VectorFunction& gradient);

  static GridSample<DataVector> fromFiniteDifferences(std::shared_ptr<Grid>&, ScalarFunction& func,
                                                       double h);

  AsInfo evaluate(AsQuadInput& input) override;

};

}  // namespace base
}  // namespace sgpp

#endif
