#include <iostream>

#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/base/tools/dimension/PcaFuncReducer.hpp>
#include <sgpp/base/tools/dimension/PcaReducer.hpp>
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/tools/dimension/AsMcReducer.hpp"
#include "sgpp/base/tools/dimension/AsQuadReducer.hpp"

double f(const sgpp::base::DataVector& v) { return v[0]; }

void f_gradient(const sgpp::base::DataVector& x, sgpp::base::DataVector& out) {
  out[0] = 1;
  out[1] = 0;
}

int main(int argc, char* argv[]) {
  auto func = sgpp::base::WrapperScalarFunction(2, f);
  size_t dim = 2;

  // Create the grid object.
  std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
  grid->getGenerator().regular(4);

  // Sample the gradient function at the grid points
  sgpp::base::SGridSample sample(grid, func);
  sample.hierarchise();

  auto solver = std::make_shared<sgpp::base::PcaCovarianceSolver>();
  auto reducerPca = sgpp::base::PcaFuncReducer(solver, std::mt19937_64::default_seed, 1000, 0.25, 100);
  sgpp::base::PcaFuncInput pcaInput{func, sample};
  sgpp::base::PcaFuncInfo info = reducerPca.evaluate(pcaInput);

  // Print out all the information that the reducer has gathered
  for (size_t d = 0; d < dim; ++d) {
    std::cout << "dimension: " + std::to_string(d) << std::endl;

    sgpp::base::DataVector e(dim);
    info.basis.getColumn(d, e);
    std::cout << "eigen vector: " + e.toString() << std::endl;
  }

  // Create the cutter used to cut off some dimensions
  auto cutter = sgpp::base::PcaFuncFixedCutter(1);

  // Use the cutter to remove unwanted dimensions
  sgpp::base::PcaFuncResult result = cutter.cut(pcaInput, info);
  std::cout << "transformation matrix used to reduce the function: " +
                   result.getProjection().getTransformationMatrix().toString()
            << std::endl;

  std::cout << "reduced dimensions: " << result.getProjection().getNewDimensions() << std::endl;

  sgpp::base::L2McRule l2Rule(std::mt19937_64::default_seed, 10000);
  double l2Error = std::sqrt(result.calculateAbsoluteError(l2Rule));
  std::cout << "L2 error: " << l2Error << std::endl;

  return 0;
}
