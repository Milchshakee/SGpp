#include "AsQuadReducer.hpp"
#include "Tools.hpp"

sgpp::base::AsQuadReducer::AsQuadReducer(std::shared_ptr<Grid>& grid,
                                         std::shared_ptr<OperationQuadrature>& quad)
    : grid(grid), quad(quad) {}

sgpp::base::AsInfo sgpp::base::AsQuadReducer::evaluate(PointSample<DataVector>& input) {

  std::vector<PointSample<double>> samples(input.getDimensions());
    for (size_t i = 0; i < input.getDimensions(); ++i) {
    std::function<double(const DataVector&)> f = [i](const DataVector& v) { return v[i]; };
      PointSample<double> grad(input.getValues(), f);
    samples[i] = grad;
      }

    sgpp::base::DataMatrix m(input.getDimensions(), input.getDimensions());
  for (size_t i = 0; i < input.getDimensions(); ++i) {
        for (size_t j = 0; j < input.getDimensions(); ++j) {
      size_t counter = 0;
          std::function<double(const DataVector&)> f = [&counter, i, j, samples](const DataVector& v) {
            double val = samples[i].getValues()[counter] * samples[j].getValues()[counter];
        counter++;
        return val;
          };
      GridSample<double> alphas(grid, f);
          DataVector v(alphas.getValues());
          double val = quad->doQuadrature(v);
          m.set(i, j, val);
          m.set(j, i, val);
    }
    }

  AsInfo i;
  i.eigenVectors = sgpp::base::DataMatrix(input.getDimensions(), input.getDimensions());
    i.eigenValues = sgpp::base::DataVector(input.getDimensions());
    Tools::svd(m, i.eigenVectors, i.eigenValues);
    return i;

}