// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DIMREDUCTION_HPP
#define DIMREDUCTION_HPP
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/tools/Sample.hpp>

namespace sgpp {
namespace base {

template <class INPUT, class INFO, class OUTPUT>
class Cutter {
 public:
  virtual OUTPUT cut(const INPUT& input, const INFO& info) = 0;
};

class ErrorRule {
public:
  virtual double calculateBasisFunction(double alpha, const HashGridPoint& point) = 0;
  virtual double calculateRelativeError(ScalarFunction& f, VectorFunction& t,
                                        ScalarFunction& r) = 0;
  virtual double calculateAbsoluteError(ScalarFunction& f, VectorFunction& t,
                                        ScalarFunction& r) = 0;
};

  class VarianceMcL2Rule : public ErrorRule {
 public:
    VarianceMcL2Rule(uint64_t seed, size_t samples);
  double calculateRelativeError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r);
    double calculateAbsoluteError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r);

  private:
  uint64_t seed;
   size_t samples;
  };

    template <class INPUT, class INFO, class OUTPUT>
  class FixedCutter : public Cutter<INPUT, INFO, OUTPUT> {
   public:
    FixedCutter(ErrorRule& r, size_t n)
      : r(r),
        n(n) {
    }

   protected:
    ErrorRule& r;
    size_t n;
  };

template <class T>
class Result {
 public:
  double calculateRelativeError(ErrorRule& c)
  {
    return c.calculateRelativeError(getOriginalFunction(), getTransformationFunction(),
                            getReducedFunction());
  }

    double calculateAbsoluteError(ErrorRule& c) {
    return c.calculateRelativeError(getOriginalFunction(), getTransformationFunction(),
                            getReducedFunction());
  }

  virtual ScalarFunction& getOriginalFunction() = 0;
  virtual ScalarFunction& getReducedFunction() = 0;
  virtual VectorFunction& getTransformationFunction() = 0;
  virtual T& getReducedOutput() = 0;
};

  template <class INPUT, class INFO, class OUTPUT>
  class ErrorRuleCutter : public Cutter<INPUT, INFO, OUTPUT> {
 public:
  ErrorRuleCutter(ErrorRule& r, double maxError) : r(r), maxError(maxError) {};

 protected:
  ErrorRule& r;
  double maxError;
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
      ProjectionFunction() = default;
      ProjectionFunction(InputProjection& p);
      ~ProjectionFunction() = default;
      void eval(const DataVector& in, DataVector& out) override;
      void clone(std::unique_ptr<VectorFunction>& clone) const override;

     private:
      InputProjection* p;
    };

    ProjectionFunction func;
  };
}  // namespace base
}  // namespace sgpp

#endif
