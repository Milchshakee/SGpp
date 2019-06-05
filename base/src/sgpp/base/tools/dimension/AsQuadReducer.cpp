#include "AsQuadReducer.hpp"
#include "Tools.hpp"

sgpp::base::AsQuadReducer::AsQuadReducer(std::shared_ptr<Grid>& grid,
                                         std::shared_ptr<OperationQuadrature>& quad)
    : grid(grid), quad(quad) {}

void sgpp::base::AsQuadReducer::evaluate(Sample<DataVector>& input, AsInfo& out) {

  std::vector<Sample<double>> samples(input.getDimensions());
    for (size_t i = 0; i < input.getDimensions(); ++i) {
    std::function<double(DataVector&)> f = [i](DataVector& v) { return v[i]; };
    Sample<double> grad = Sample<double>(input.getValues(), f);
    samples[i] = grad;
      }

    sgpp::base::DataMatrix m(input.getDimensions(), input.getDimensions());
  for (size_t i = 0; i < input.getDimensions(); ++i) {
        for (size_t j = 0; j < input.getDimensions(); ++j) {
      size_t counter = 0;
          std::function<double(DataVector&)> f = [&counter, i, j, samples](DataVector& v) {
            double val = samples[i].getValues()[counter] * samples[j].getValues()[counter];
        counter++;
        return val;
          };
      Sample<double> alphas = Sampler::sampleGrid<double>(*grid, f);
          DataVector v(alphas.getValues());
          double val = quad->doQuadrature(v);
          m.set(i, j, val);
          m.set(j, i, val);
    }
    }

  
  out.eigenVectors = sgpp::base::DataMatrix(input.getDimensions(), input.getDimensions());
    out.eigenValues = sgpp::base::DataVector(input.getDimensions());
    Tools::svd(m, out.eigenVectors, out.eigenValues);

}

sgpp::base::AsResult sgpp::base::AsQuadReducer::reduce(Sample<DataVector>& input, size_t c,
  const AsInfo& info) {
}
