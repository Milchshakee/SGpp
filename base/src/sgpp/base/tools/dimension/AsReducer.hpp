// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASREDUCER_HPP
#define ASREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include "DimReduction.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {

struct AsInfo {
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
};

class AsReducedFunction : public sgpp::optimization::ScalarFunction {
 public:
  AsReducedFunction(std::unique_ptr<sgpp::optimization::ScalarFunction>&& function,
                    DataMatrix transformation);
  ~AsReducedFunction() override = default;

  double eval(const base::DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;

 private:
  std::unique_ptr<sgpp::optimization::ScalarFunction> function;
  DataMatrix transformation;
};

struct AsResult {
  AsResult(DataMatrix& m, size_t n);

  DataMatrix transformation;

  std::unique_ptr<AsReducedFunction> apply(sgpp::optimization::ScalarFunction& input);
};

  template <class INPUT>
class AsEigenValueCutter : public Cutter<INPUT, AsInfo, AsResult> {
 public:
    AsEigenValueCutter(double minValue) : minEigenValue(minValue) {}

  void cut(const INPUT& input, const AsInfo& info) override
  {
    for (size_t d = 0; d < info.eigenValues.size(); ++d) {
      if (info.eigenValues[d] < minEigenValue) {
        return d + 1;
      }
    }
    return info.eigenValues.size();
  };

 private:
  double minEigenValue;
};

template <class I>
class AsReducer : public sgpp::base::Reducer<I, AsInfo, AsResult> {
};

}  // namespace base
}  // namespace sgpp

#endif
