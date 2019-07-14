#include "AsQuadReducer.hpp"
#include "EigenHelper.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"


sgpp::base::AsQuadFixedCutter::AsQuadFixedCutter(size_t n) : n(n) {
}

sgpp::base::AsResult sgpp::base::AsQuadFixedCutter::cut(const GridSample<DataVector>& input,
                                                        const AsInfo& info) {
  return AsResult(info.eigenVectors, n);
}

sgpp::base::AsQuadEigenValueCutter::AsQuadEigenValueCutter(double minEigenValue) : AsEigenValueCutter<sgpp::base::GridSample<sgpp::base::DataVector>>(minEigenValue) {
}

sgpp::base::AsInfo sgpp::base::AsQuadReducer::evaluate(GridSample<DataVector>& input) {

  std::vector<PointSample<double>> samples(input.getDimensions());
    for (size_t i = 0; i < input.getDimensions(); ++i) {
    std::function<double(const DataVector&)> f = [i](const DataVector& v) { return v[i]; };
      PointSample<double> grad(input.getValues(), f);
    samples[i] = grad;
      }
    std::unique_ptr<sgpp::base::OperationQuadrature> op(
        sgpp::op_factory::createOperationQuadrature(const_cast<Grid&>(input.getGrid())));
    sgpp::base::DataMatrix m(input.getDimensions(), input.getDimensions());
  for (size_t i = 0; i < input.getDimensions(); ++i) {
        for (size_t j = 0; j < input.getDimensions(); ++j) {
      size_t counter = 0;
          std::function<double(const DataVector&)> f = [&counter, i, j, samples](const DataVector& v) {
            double val = samples[i].getValues()[counter] * samples[j].getValues()[counter];
        counter++;
        return val;
          };
      PointSample<double> alphas(input.getKeys(), f);
          DataVector v(alphas.getValues());
          double val = op->doQuadrature(v);
          m.set(i, j, val);
          m.set(j, i, val);
    }
    }

  AsInfo i;
  i.eigenVectors = sgpp::base::DataMatrix(input.getDimensions(), input.getDimensions());
    i.eigenValues = sgpp::base::DataVector(input.getDimensions());
    EigenHelper::svd(EigenHelper::toEigen(m), i.eigenVectors, i.eigenValues);
    return i;

}
