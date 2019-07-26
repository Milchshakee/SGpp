#include "AsQuadReducer.hpp"
#include "EigenHelper.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"


sgpp::base::AsQuadResult::AsQuadResult(const SGridSample& input, const DataMatrix& m, size_t n) : AsResult<sgpp::base::SGridSample>(m, n) {

}

sgpp::optimization::ScalarFunction& sgpp::base::AsQuadResult::getReducedFunction() {

}

sgpp::base::SGridSample& sgpp::base::AsQuadResult::getReducedOutput() { return reduced; }

sgpp::base::AsQuadFixedCutter::AsQuadFixedCutter(size_t n) : n(n) {
}

sgpp::base::AsQuadResult sgpp::base::AsQuadFixedCutter::cut(const AsQuadInput& input,
                                                        const AsInfo& info) {
  return AsQuadResult(input.functionSample, info.eigenVectors, n);
}

sgpp::base::AsInfo sgpp::base::AsQuadReducer::evaluate(AsQuadInput& input) {

  std::vector<PointSample<double>> samples(input.gradientSample.getDimensions());
    for (size_t i = 0; i < input.gradientSample.getDimensions(); ++i) {
    std::function<double(const DataVector&)> f = [i](const DataVector& v) { return v[i]; };
      PointSample<double> grad(input.gradientSample.getValues(), f);
    samples[i] = grad;
      }
    std::unique_ptr<sgpp::base::OperationQuadrature> op(
        sgpp::op_factory::createOperationQuadrature(const_cast<Grid&>(input.gradientSample.getGrid())));
    sgpp::base::DataMatrix m(input.gradientSample.getDimensions(), input.gradientSample.getDimensions());
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
