// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DIMREDUCTION_HPP
#define DIMREDUCTION_HPP
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/Sample.hpp>

namespace sgpp {
namespace base {

template <class INPUT, class INFO, class OUTPUT>
class Cutter {
 public:
  virtual OUTPUT cut(const INPUT& input, const INFO& info) = 0;
};

class Coverage {
public:
  virtual double calculateError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r) = 0;
};

template <class T>
class Result {
 public:
  double calcMcL2Error(ScalarFunction& func, size_t paths,
                       uint64_t seed = std::mt19937_64::default_seed)
  {
    std::mt19937_64 rand(seed);
    std::uniform_real_distribution<double> dist(0, 1);
    size_t funcDimensions = func.getNumberOfParameters();
    size_t newDimensions = getReducedFunction().getNumberOfParameters();

    sgpp::base::DataVector point(funcDimensions);
    double res = 0;

    for (size_t i = 0; i < paths; i++) {
      for (size_t d = 0; d < funcDimensions; d++) {
        point[d] = dist(rand);
      }
      double val = func.eval(point);
      DataVector out(newDimensions);
      getTransformationFunction().eval(point, out);
      res += pow(val - getReducedFunction().eval(out), 2);
    }

    return sqrt(res / static_cast<double>(paths));
  }

  double getError(Coverage& c)
  {
    return c.calculateError(getOriginalFunction(), getTransformationFunction(),
                            getReducedFunction());
  }

  virtual double getCoveredVariance() = 0;
  virtual ScalarFunction& getOriginalFunction() = 0;
  virtual ScalarFunction& getReducedFunction() = 0;
  virtual VectorFunction& getTransformationFunction() = 0;
  virtual T& getReducedOutput() = 0;
};

template <class INPUT, class INFO, class OUTPUT>
class Reducer {
 public:
  Reducer() = default;
  virtual ~Reducer() = default;

  OUTPUT evaluateAndCut(INPUT& input, Cutter<INPUT, INFO, OUTPUT>& cutter) {
    INFO i = evaluateFunction(input);
    OUTPUT result = cutter.cut(input, i);
    return result;
  }

  virtual INFO evaluate(INPUT& input) = 0;
};

  class InputProjection
  {
  public:
    InputProjection(const DataMatrix& basis, size_t n, const DataVector& mean);

    void inverse(const DataVector& in, DataVector& out);
    VectorFunction& getFunction();

  private:
    DataMatrix basis;
    DataMatrix cutBasis;
    size_t oldDimensions;
    size_t newDimensions;
    DataVector mean;
    DataVector posRange;
    DataVector negRange;
    DataVector start;

    void calculateRanges();

    class ProjectionFunction : public VectorFunction {
     public:
      ProjectionFunction(InputProjection& p);
      ~ProjectionFunction() = default;
      void eval(const DataVector& in, DataVector& out) override;
      void clone(std::unique_ptr<VectorFunction>& clone) const override;

     private:
      InputProjection& p;
    };

    ProjectionFunction func;
  };
}  // namespace base
}  // namespace sgpp

#endif
