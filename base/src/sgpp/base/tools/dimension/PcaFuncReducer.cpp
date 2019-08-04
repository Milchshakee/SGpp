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
    : mean(mean), newDimensions(n), oldDimensions(m.getNrows()) {
  basis = m;
  this->mInv = m;
  this->mInv.resizeRowsCols(m.getNrows(), n);
  DataMatrix copy = m;
  copy.transpose();

  calculateRanges();

  std::shared_ptr<Grid> newGrid(const_cast<Grid&>(input.getGrid()).createGridOfEquivalentType(n));
  newGrid->getGenerator().regular(const_cast<Grid&>(input.getGrid()).getStorage().getMaxLevel());

  EvalFunction e(input);
  size_t dim = input.getDimensions();
  std::function<double(const DataVector&)> func = [this, &e, dim, n](const DataVector& v) {
    DataVector newV;
    transformFrom(v, newV);
    std::cout << "e: " + newV.toString() << std::endl;
    return e.eval(newV);
  };

  SGridSample newSample(newGrid, func);
  newSample.hierarchise();
  reduced = newSample;
  evalFunc = EvalFunction(reduced);
}

void sgpp::base::PcaFuncResult::calculateRanges() {
  size_t dim = newDimensions;
  posRange = DataVector(dim, 2);
  negRange = DataVector(dim, -2);
  for (size_t d = 0; d < dim; ++d) {
    double m = mean[d];
    for (size_t i = 0; i < oldDimensions; ++i) {
      DataVector col(oldDimensions);
      basis.getColumn(i, col);
      double v = col[d];
      double posX = (1 - m) / v;
      double negX = (0 - m) / v;
      if (posX < 0) {
        std::swap(posX, negX);
      }
      if (posX < posRange[d]) {
        posRange[d] = posX;
      }
      if (negX > negRange[d]) {
        negRange[d] = negX;
      }
    }
  }
}

sgpp::base::ScalarFunction& sgpp::base::PcaFuncResult::getReducedFunction() {}

sgpp::base::VectorFunction& sgpp::base::PcaFuncResult::getTransformationFunction() {}

sgpp::base::SGridSample& sgpp::base::PcaFuncResult::getReducedOutput() {}

void sgpp::base::PcaFuncResult::transformFrom(const DataVector& in, DataVector& out) {
  DataVector start = mean;
  for (size_t d = 0; d < newDimensions; ++d) {
    DataVector add(oldDimensions);
    basis.getColumn(d, add);
    add.mult(negRange[d]);
    start.add(add);
  }

  for (size_t d = 0; d < newDimensions; ++d) {
    double scale = posRange[d] - negRange[d];
    DataVector add(oldDimensions);
    basis.getColumn(d, add);
    add.mult(scale * in[d]);
    start.add(add);
  }

  out = start;
}

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
  double min = 0;
  for (size_t i = 0; i < input.getSize(); ++i) {
    double val = input.getValues()[i];
    if (val < min) {
      val = min;
    }
  }

  SGridSample copy = input;
  if (min < 0) {
    for (size_t i = 0; i < input.getSize(); ++i) {
      copy.getValues()[i] -= min;
    }
  }

  copy.hierarchise();
  double quad = copy.quadrature();
  double scale = 1.0 / quad;
  for (size_t i = 0; i < input.getSize(); ++i) {
    input.getValues()[i] *= scale;
  }
  input.hierarchise();
}
