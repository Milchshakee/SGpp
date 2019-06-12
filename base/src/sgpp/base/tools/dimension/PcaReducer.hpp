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
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
};

  struct PcaResult
  {
  PcaResult(const DataMatrix& m, size_t n);

  DataMatrix transformation;

  FixedDistribution apply(const VectorDistribution& input);
  };

  class PcaCutter : public Cutter<VectorDistribution, PcaInfo, PcaResult>
  {
  };

    class PcaVarianceCutter : public PcaCutter {
   public:
      PcaVarianceCutter(double variancePercentage);

    PcaResult cut(const VectorDistribution& input, const PcaInfo& info) override;

    private:
    double variancePercentage;
  };

  class PcaFixedCutter : public PcaCutter
  {
  public:
    PcaFixedCutter(size_t n);

   PcaResult cut(const VectorDistribution& input, const PcaInfo& info) override;

  private:
   size_t n;
  };

class PcaReducer : public Reducer<VectorDistribution, PcaInfo, PcaResult> {
 public:
  PcaReducer(size_t iterations, uint64_t seed);

  PcaInfo evaluate(VectorDistribution& input) override;

 private:
  size_t iterations;
  uint64_t seed;
};

}  // namespace base
}  // namespace sgpp

#endif
