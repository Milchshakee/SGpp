// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ANOVAREDUCER_HPP
#define ANOVAREDUCER_HPP

#include <sgpp/base/tools/dimension/DimReduction.hpp>

/**
 * Vector that holds levels for every dimension
 */
typedef std::vector<bool> DimensionVector;

struct AnovaInformation
{
  struct Component
  {
    size_t order;
    double variance;
    DimensionVector fixedDimensions;
  };

  std::vector<double> variances;
  std::vector<Component> components;
};

typedef std::vector<DimensionVector> AnovaCutoff;

class ANOVAReducer : public sgpp::base::FunctionReducer<AnovaInformation, AnovaCutoff> {
 public:
  enum class GridType { FULL, SPARSE };

  class VarianceCutoff : public sgpp::base::CutoffCriterion<AnovaInformation, AnovaCutoff> {
   public:
    VarianceCutoff(double maxVariance);

    AnovaCutoff evaluate(const AnovaInformation& info) override;

   private:
    double maxVariance;
  };

    class OrderCutoff : public sgpp::base::CutoffCriterion<AnovaInformation, AnovaCutoff> {
   public:
    OrderCutoff(double maxVariance);

    AnovaCutoff evaluate(const AnovaInformation& info) override;

   private:
    double maxVariance;
  };

  ANOVAReducer(size_t gridLevel, GridType gridType);

  void evaluateFunction(sgpp::optimization::ScalarFunction& input, AnovaInformation& out) override;
  std::unique_ptr<sgpp::optimization::ScalarFunction> reduceFunction(
      sgpp::optimization::ScalarFunction& input, const AnovaCutoff& c, const AnovaInformation& info) override;

 private:
  size_t gridLevel;
  GridType gridType;
};

#endif
