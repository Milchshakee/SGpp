#include <sgpp/base/function/scalar/EvalFunction.hpp>
#include <sgpp/base/tools/dimension/AsQuadReducer.hpp>


sgpp::base::GridSample<sgpp::base::DataVector> sgpp::base::AsQuadReducer::fromGradientFunction(
  std::shared_ptr<Grid>& grid, VectorFunction& gradient) {
  sgpp::base::GridSample<sgpp::base::DataVector> gradientSamples =
      sgpp::base::SampleHelper::sampleGrid(grid, gradient);
  return gradientSamples;
}

sgpp::base::GridSample<sgpp::base::DataVector> sgpp::base::AsQuadReducer::fromFiniteDifferences(
    std::shared_ptr<Grid>& grid, ScalarFunction& func, double h)
{
  std::vector<DataVector> samples(grid->getSize(), DataVector(grid->getDimension()));
  GridDistribution dist(*grid);
  DataVector working;
  for (size_t i = 0; i < grid->getSize(); ++i) {
    const DataVector& vec = dist.getVectors()[i];
    working = vec;
    for (size_t d = 0; d < grid->getDimension(); ++d) {
      double hOffset = vec[d] + h;
      if (hOffset > 1.0) {
        hOffset -= 2 * h;
      }
      working[d] = hOffset;
      double val = (func.eval(vec) - func.eval(working)) / h;
      samples[i][d] = val;
      working[d] = vec[d];
    }
  }
  return GridSample<DataVector>(grid, samples);
}

sgpp::base::AsQuadResult::AsQuadResult(const AsQuadInput& input, const DataMatrix& m, size_t n,
                                       const DataVector& mean)
    : projection(m, n, mean), originalFunction(&input.originalFunction) {
  std::shared_ptr<Grid> newGrid(const_cast<Grid&>(input.functionSample.getGrid()).createGridOfEquivalentType(n));
  newGrid->getGenerator().regular(
      const_cast<Grid&>(input.functionSample.getGrid()).getStorage().getMaxLevel());

  std::function<double(const DataVector&)> func = [this,&input](const DataVector& v) {
    DataVector newV;
    projection.inverse(v, newV);
    return input.originalFunction.eval(newV);
  };

  SGridSample newSample(newGrid, func);
  newSample.hierarchise();
  reduced = newSample;
  evalFunc = EvalFunction(reduced);
}

sgpp::base::ScalarFunction& sgpp::base::AsQuadResult::getReducedFunctionSurrogate() {
  return evalFunc;
}

sgpp::base::SGridSample& sgpp::base::AsQuadResult::getReducedOutput() { return reduced; }


sgpp::base::ScalarFunction& sgpp::base::AsQuadResult::getOriginalFunction() { return *originalFunction; }

sgpp::base::VectorFunction& sgpp::base::AsQuadResult::getTransformationFunction() {
  return projection.getFunction();
}


sgpp::base::InputProjection& sgpp::base::AsQuadResult::getProjection() { return projection; }

sgpp::base::AsQuadFixedCutter::AsQuadFixedCutter(size_t n, const DataVector& mean)
    : FixedCutter<sgpp::base::AsQuadInput, sgpp::base::AsInfo, sgpp::base::AsQuadResult>(n), mean(mean) {}

sgpp::base::AsQuadResult sgpp::base::AsQuadFixedCutter::cut(const AsQuadInput& input,
                                                            const AsInfo& info) {
  return AsQuadResult(input, info.eigenVectors, n, mean);
}


sgpp::base::AsQuadErrorRuleCutter::AsQuadErrorRuleCutter(ErrorRule& rule, double minValue,
  const DataVector& mean) : ErrorRuleCutter<sgpp::base::AsQuadInput, sgpp::base::AsInfo, sgpp::base::AsQuadResult>(rule, minValue), mean(mean) {
}

sgpp::base::AsQuadResult sgpp::base::AsQuadErrorRuleCutter::cut(const AsQuadInput& input,
                                                                const AsInfo& info) {
  size_t dim = input.originalFunction.getNumberOfParameters();
  AsQuadResult last(input, info.eigenVectors, dim,
                  mean);
  for (size_t d = 1; d < dim; ++d) {
    AsQuadResult result(input, info.eigenVectors, dim - d, mean);
    double var = 1 - result.calculateRelativeError(r);
    if (var > minVariance) {
      last = result;
    } else {
      return last;
    }
  }
  return last;
}

sgpp::base::AsInfo sgpp::base::AsQuadReducer::evaluate(AsQuadInput& input) {
  std::vector<PointSample<double>> samples(input.gradientSample.getDimensions());
  for (size_t i = 0; i < input.gradientSample.getDimensions(); ++i) {
    std::function<double(const DataVector&)> f = [i](const DataVector& v) { return v[i]; };
    PointSample<double> grad(input.gradientSample.getValues(), f);
    samples[i] = grad;
  }
  std::unique_ptr<sgpp::base::OperationQuadrature> op(sgpp::op_factory::createOperationQuadrature(
      const_cast<Grid&>(input.gradientSample.getGrid())));
  sgpp::base::DataMatrix m(input.gradientSample.getDimensions(),
                           input.gradientSample.getDimensions());
  for (size_t i = 0; i < input.gradientSample.getDimensions(); ++i) {
    for (size_t j = 0; j < input.gradientSample.getDimensions(); ++j) {
      size_t counter = 0;
      std::function<double(const DataVector&)> f = [&counter, i, j, samples](const DataVector& v) {
        double val = samples[i].getValues()[counter] * samples[j].getValues()[counter];
        counter++;
        return val;
      };
      PointSample<double> alphas(input.gradientSample.getKeys(), f);
      DataVector v(alphas.getValues());
      double val = op->doQuadrature(v);
      m.set(i, j, val);
      m.set(j, i, val);
    }
  }

  AsInfo i;
  i.eigenVectors = sgpp::base::DataMatrix(input.gradientSample.getDimensions(),
                                          input.gradientSample.getDimensions());
  i.eigenValues = sgpp::base::DataVector(input.gradientSample.getDimensions());
  EigenHelper::svd(EigenHelper::toEigen(m), i.eigenVectors, i.eigenValues);
  return i;
}
