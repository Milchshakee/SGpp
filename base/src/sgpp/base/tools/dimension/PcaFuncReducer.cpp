#include <sgpp/base/function/scalar/EvalFunction.hpp>
#include <sgpp/base/tools/dimension/PcaFuncReducer.hpp>
#include <sgpp/base/tools/dist/RandomPdfDistribution.hpp>

sgpp::base::PcaFuncFixedCutter::PcaFuncFixedCutter(size_t n)
    : FixedCutter<sgpp::base::PcaFuncInput, sgpp::base::PcaFuncInfo, sgpp::base::PcaFuncResult>(n) {
}

sgpp::base::PcaFuncResult sgpp::base::PcaFuncFixedCutter::cut(const PcaFuncInput& input,
                                                              const PcaFuncInfo& info) {
  return PcaFuncResult(input, info.basis, n, info.mean);
}

sgpp::base::PcaFuncResult sgpp::base::PcaFuncErrorRuleCutter::cut(const PcaFuncInput& input,
                                                                  const PcaFuncInfo& info) {
  PcaFuncResult last(input, info.basis, input.sample.getDimensions(), info.mean);
  for (size_t d = 1; d < input.sample.getDimensions(); ++d) {
    PcaFuncResult result(input, info.basis, input.sample.getDimensions() - d, info.mean);
    if (result.calculateRelativeError(r) > maxError) {
      return last;
    } else {
      last = result;
    }
  }
  return last;
}

sgpp::base::PcaFuncResult::PcaFuncResult(const PcaFuncInput& input, const DataMatrix& m, size_t n,
                                         const DataVector& mean)
    : projection(m, n, mean), originalFunction(&input.function) {
  std::shared_ptr<Grid> newGrid(
      const_cast<Grid&>(input.sample.getGrid()).createGridOfEquivalentType(n));
  newGrid->getGenerator().regular(
      const_cast<Grid&>(input.sample.getGrid()).getStorage().getMaxLevel());

  std::function<double(const DataVector&)> func = [this, &input](const DataVector& v) {
    DataVector newV;
    projection.inverse(v, newV);
    return input.function.eval(newV);
  };

  SGridSample newSample(newGrid, func);
  newSample.hierarchise();
  reducedSurrogate = newSample;
  reducedSurrogateFunction = EvalFunction(reducedSurrogate);
}


sgpp::base::ScalarFunction& sgpp::base::PcaFuncResult::getOriginalFunction() {
  return *originalFunction;
}

sgpp::base::ScalarFunction& sgpp::base::PcaFuncResult::getReducedFunctionSurrogate() { return reducedSurrogateFunction; }

sgpp::base::VectorFunction& sgpp::base::PcaFuncResult::getTransformationFunction() {
  return projection.getFunction();
}

sgpp::base::SGridSample& sgpp::base::PcaFuncResult::getReducedOutput() { return reducedSurrogate; }


sgpp::base::InputProjection& sgpp::base::PcaFuncResult::getProjection() { return projection; }

sgpp::base::PcaFuncReducer::PcaFuncReducer(std::shared_ptr<PcaSolver> solver, uint64_t seed,
                                           size_t samples, double stepSize, size_t iterations)
    : solver(solver), seed(seed), samples(samples), stepSize(stepSize), iterations(iterations) {}

sgpp::base::PcaFuncInfo sgpp::base::PcaFuncReducer::evaluate(PcaFuncInput& input) {
  if (!input.sample.isHierarchised()) {
    throw std::invalid_argument("Input is not hierarchised");
  }

  SGridSample density = input.sample;
  toDensity(density);
  EvalFunction f(density);
  //ScalarFunction& f = input.function;
  RandomPdfDistribution d(samples, input.sample.getDimensions(), seed, f, iterations, stepSize);
  auto reducer = sgpp::base::PcaReducer(solver);
  sgpp::base::PcaInfo info = reducer.evaluate(d);

  size_t dim = input.sample.getDimensions();
  PcaFuncInfo toReturn;
  toReturn.mean = info.mean;
  toReturn.basis = DataMatrix(dim, dim);

  for (size_t d = 0; d < dim; ++d) {
    sgpp::base::DataVector pa(input.sample.getDimensions());
    info.basis.getColumn(d, pa);
    toReturn.basis.setColumn(dim - d - 1, pa);
  }
  return toReturn;
}

void sgpp::base::PcaFuncReducer::toDensity(SGridSample& input) {
  input.dehierarchise();

  double min = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < input.getSize(); ++i) {
    double val = input.getValues()[i];
    if (val < min) {
      val = min;
    }
  }

  for (size_t i = 0; i < input.getSize(); ++i) {
    input.getValues()[i] -= min;
  }

  input.hierarchise();
  double quad = input.quadrature();
  if (quad == 0.0) {
    for (size_t i = 0; i < input.getSize(); ++i) {
      input.getValues()[i] = 1;
    }
  } else {
    double scale = 1.0 / quad;
    for (size_t i = 0; i < input.getSize(); ++i) {
      input.getValues()[i] *= scale;
    }
  }
}
