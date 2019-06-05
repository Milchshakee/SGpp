// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ACTIVESUBSPACEREDUCER_HPP
#define ACTIVESUBSPACEREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "AsReducer.hpp"
#include "sgpp/base/operation/hash/OperationQuadrature.hpp"

namespace sgpp {
namespace base {

class AsQuadReducer : public AsReducer<Sample<DataVector>> {
 public:

  AsQuadReducer(std::shared_ptr<Grid>& grid, std::shared_ptr<OperationQuadrature>& quad);

  void evaluate(Sample<DataVector>& input, AsInfo& out) override;
  AsResult reduce(Sample<DataVector>& input, size_t c,
                              const AsInfo& info) override;
private:
  std::shared_ptr<Grid> grid;
 std::shared_ptr<OperationQuadrature> quad;
};

}  // namespace base
}  // namespace sgpp

#endif
