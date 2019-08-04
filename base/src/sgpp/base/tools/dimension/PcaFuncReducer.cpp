#include "PcaFuncReducer.hpp"
#include <sgpp/base/function/scalar/EvalFunction.hpp>
#include <sgpp/base/tools/dist/RandomPdfDistribution.hpp>
#include "EigenHelper.hpp"

sgpp::base::PcaFuncFixedCutter::PcaFuncFixedCutter(size_t n) : n(n) {}

sgpp::base::PcaFuncResult sgpp::base::PcaFuncFixedCutter::cut(const SGridSample& input,
                                                              const PcaFuncInfo& info) {
  return PcaFuncResult(input, info.basis, n, info.mean);
}

sgpp::base::PcaFuncVarianceCutter::PcaFuncVarianceCutter(double minVarianceShare) {}

sgpp::base::PcaFuncResult sgpp::base::PcaFuncVarianceCutter::cut(const SGridSample& input,
                                                                 const PcaFuncInfo& info) {}

sgpp::base::PcaFuncResult::PcaFuncResult(const SGridSample& input, const DataMatrix& m, size_t n,
                                         const DataVector& mean)
    : projection(m, n, mean) {
  std::shared_ptr<Grid> newGrid(const_cast<Grid&>(input.getGrid()).createGridOfEquivalentType(n));
  newGrid->getGenerator().regular(const_cast<Grid&>(input.getGrid()).getStorage().getMaxLevel());

  EvalFunction e(input);
  size_t dim = input.getDimensions();
  std::function<double(const DataVector&)> func = [this, &e, dim, n](const DataVector& v) {
    DataVector newV;
    projection.inverse(v, newV);
    std::cout << "e: " + newV.toString() << std::endl;
    return e.eval(newV);
  };

  SGridSample newSample(newGrid, func);
  newSample.hierarchise();
  reduced = newSample;
  evalFunc = EvalFunction(reduced);
}

sgpp::base::ScalarFunction& sgpp::base::PcaFuncResult::getReducedFunction()
{ return evalFunc; }

sgpp::base::VectorFunction& sgpp::base::PcaFuncResult::getTransformationFunction()
{
  return projection.getFunction();
}

sgpp::base::SGridSample& sgpp::base::PcaFuncResult::getReducedOutput()
{ return reduced; }

sgpp::base::PcaFuncReducer::PcaFuncReducer(std::shared_ptr<PcaSolver> solver, uint64_t seed,
                                           size_t samples, double stepSize, size_t iterations)
    : solver(solver), seed(seed), samples(samples), stepSize(stepSize), iterations(iterations) {}

sgpp::base::PcaFuncInfo sgpp::base::PcaFuncReducer::evaluate(SGridSample& input) {
  SGridSample density = input;
  toDensity(density);
  EvalFunction f(density);
  RandomPdfDistribution d(samples, input.getDimensions(), seed, f, iterations, stepSize);
  auto reducer = sgpp::base::PcaReducer(solver);
  sgpp::base::PcaInfo info = reducer.evaluate(d);

  size_t a = info.activeComponentsCount;
  size_t dim = input.getDimensions();
  PcaFuncInfo toReturn;
  toReturn.mean = info.mean;
  toReturn.varianceShares = DataVector(dim);
  toReturn.basis = DataMatrix(dim, dim);
  toReturn.eigenValues = DataVector(dim);

  for (size_t d = 0; d < dim; ++d) {
    toReturn.varianceShares[d] = d >= a ? 0 : info.varianceShares[a - d - 1];
    sgpp::base::DataVector pa(input.getDimensions());
    info.basis.getColumn(d, pa);
    toReturn.basis.setColumn(dim - d - 1, pa);
    toReturn.eigenValues.set(dim - d - 1, info.eigenValues.get(d));
  }
  return toReturn;
}

void sgpp::base::PcaFuncReducer::toDensity(SGridSample& input) {
  double min = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < input.getSize(); ++i) {
    double val = input.getValues()[i];
    if (val < min) {
      val = min;
    }
  }

  SGridSample copy = input;
  for (size_t i = 0; i < input.getSize(); ++i) {
    copy.getValues()[i] -= min;
  }

  copy.hierarchise();
  double quad = copy.quadrature();
  double scale = 1.0 / quad;
  for (size_t i = 0; i < input.getSize(); ++i) {
    input.getValues()[i] *= scale;
  }
  input.hierarchise();
}
