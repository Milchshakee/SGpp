#include <iostream>

#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/base/tools/dimension/PcaFuncReducer.hpp>
#include <sgpp/base/tools/dimension/PcaReducer.hpp>
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/tools/dimension/AsMcReducer.hpp"
#include "sgpp/base/tools/dimension/AsQuadReducer.hpp"

double f(const sgpp::base::DataVector& v) { return v[0] + v[1]; }

int main(int argc, char* argv[]) {
  auto func = sgpp::base::WrapperScalarFunction(2, f);
  size_t dim = 2;

  // Create the grid object.
  std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
  grid->getGenerator().regular(8);

  // Sample the gradient function at the grid points
  sgpp::base::SGridSample sample(grid, func);
  sample.hierarchise();

  // Create the reducer
  auto reducer = sgpp::base::AsMcReducer();

  sgpp::base::EvalFunction sampleFunc(sample);
  sgpp::base::GridDistribution dist(*grid);
  sgpp::base::PointSample<sgpp::base::DataMatrix> m =
      sgpp::base::AsMcReducer::fromFiniteDifferences(sampleFunc, dist, 0.0001);
  sgpp::base::AsMcInput i{func, m};

  // Use the reducer to first evaluate the function using the sampled gradients
  sgpp::base::AsInfo info = reducer.evaluate(i);

  // Print out all the information that the reducer has gathered
  for (size_t d = 0; d < dim; ++d) {
    std::cout << "dimension: " + std::to_string(d) << std::endl;

    sgpp::base::DataVector e(dim);
    info.eigenVectors.getColumn(d, e);
    std::cout << "eigen vector: " + e.toString() << std::endl;

    std::cout << "eigen value: " + std::to_string(info.eigenValues[d]) << std::endl;
    std::cout << std::endl;
  }

  // Calculate the mean using PCA
  auto solver = std::make_shared<sgpp::base::PcaCovarianceSolver>();
  auto reducerPca = sgpp::base::PcaFuncReducer(solver, std::mt19937_64::default_seed, 100, 0.5, 10);
  sgpp::base::PcaFuncInput pcaInput{func, sample};
  sgpp::base::PcaFuncInfo infoPca = reducerPca.evaluate(pcaInput);

  // Create the cutter used to cut off some dimensions
  auto cutter = sgpp::base::AsMcFixedCutter(1, sgpp::base::GridType::LinearBoundary, 8, infoPca.mean);

  // Use the cutter to remove unwanted dimensions
  sgpp::base::AsMcResult result = cutter.cut(i, info);
  std::cout << "transformation matrix used to reduce the function: " +
                   result.getProjection().getTransformationMatrix().toString()
            << std::endl;

  std::cout << "reduced dimensions: " << result.getProjection().getNewDimensions() << std::endl;

  sgpp::base::L2McRule l2Rule(std::mt19937_64::default_seed, 10000);
  double l2Error = std::sqrt(result.calculateAbsoluteError(l2Rule));
  std::cout << "L2 error: " << l2Error << std::endl;

  return 0;
}
