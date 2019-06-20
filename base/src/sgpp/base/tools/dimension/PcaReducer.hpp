// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ITERATIVEPCAREDUCER_HPP
#define ITERATIVEPCAREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "AsMcReducer.hpp"

namespace sgpp {
namespace base {

  struct PcaInfo {
  sgpp::base::DataMatrix principalAxes;
  sgpp::base::DataMatrix loadings;
    sgpp::base::DataVector eigenValues;
    sgpp::base::DataVector singularValues;
  DataVector varianceShares;
  };

  struct PcaResult
  {
  PcaResult(const DataMatrix& m, size_t n, double coveredVariance);

  DataMatrix transformation;
  double coveredVariance;

  FixedDistribution apply(const VectorDistribution& input);
  };

  class PcaCutter : public Cutter<VectorDistribution, PcaInfo, PcaResult>
  {
  };

    class PcaVarianceCutter : public PcaCutter {
   public:
      PcaVarianceCutter(double minVarianceShare);

    PcaResult cut(const VectorDistribution& input, const PcaInfo& info) override;

    private:
    double minVarianceShare;
  };

  class PcaFixedCutter : public PcaCutter
  {
  public:
    PcaFixedCutter(size_t n);

   PcaResult cut(const VectorDistribution& input, const PcaInfo& info) override;

  private:
   size_t n;
  };

  class PcaSolver
  {
  public:
    virtual PcaInfo solve(DataMatrix& matrix) = 0;
  };

  class PcaCovarianceSolver : public PcaSolver {
   public:
    PcaInfo solve(DataMatrix& matrix) override;
  };

  class PcaIterativeSolver : public PcaSolver {
   public:
    PcaIterativeSolver(size_t iterations, uint64_t seed);

    PcaInfo solve(DataMatrix& matrix) override;

  private:
    size_t iterations;
   uint64_t seed;
  };

class PcaReducer : public Reducer<VectorDistribution, PcaInfo, PcaResult> {
 public:
  PcaReducer(std::shared_ptr<PcaSolver> solver);

  PcaInfo evaluate(VectorDistribution& input) override;

 private:
  std::shared_ptr<PcaSolver> solver;
};

}  // namespace base
}  // namespace sgpp

#endif
