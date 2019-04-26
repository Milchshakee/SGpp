// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ACTIVESUBSPACEREDUCER_HPP
#define ACTIVESUBSPACEREDUCER_HPP

#include <sgpp/base/tools/dimension/Reducer.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include "VectorDistribution.hpp"

class ActiveSubspaceReducer : public Reducer {
 public:
  class GradientGenerationStrategy
  {
  public:
    virtual sgpp::base::DataVector gradientAt(sgpp::base::DataVector& v);
  };

  class MatrixComputationStrategy
  {
  public:
    MatrixComputationStrategy(size_t dimensions, GradientGenerationStrategy& s);
    virtual ~MatrixComputationStrategy() = default;

    virtual std::unique_ptr<sgpp::base::DataMatrix> compute();
  protected:
    size_t dimensions;
    GradientGenerationStrategy& gradient;
  };

  class MCStrategy : MatrixComputationStrategy
  {
  public:
    MCStrategy(size_t dimensions, GradientGenerationStrategy& s, VectorDistribution& distribution,
               size_t samples);

    std::unique_ptr<sgpp::base::DataMatrix> compute() override;

  private:
    VectorDistribution distribution;
   size_t samples;
  };

  class CutoffCriterion
  {
    
  };

  ActiveSubspaceReducer(MatrixComputationStrategy matrix);

  std::unique_ptr<sgpp::optimization::VectorFunction> reduce(
      sgpp::optimization::VectorFunction& input);

private:
  MatrixComputationStrategy matrix;
};

#endif
