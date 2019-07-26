#include <iostream>

#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/tools/dimension/AsMcReducer.hpp"
#include "sgpp/base/tools/dimension/AsQuadReducer.hpp"
#include "sgpp/optimization/function/scalar/WrapperScalarFunction.hpp"
#include "sgpp/optimization/function/vector/WrapperVectorFunction.hpp"

double f(const sgpp::base::DataVector& v) { return v[0]; }

void f_gradient(const sgpp::base::DataVector& x, sgpp::base::DataVector& out) {
  out[0] = 1;
  out[1] = 0;
}


int main(int argc, char* argv[]) {
  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
  auto funcGradient = sgpp::optimization::WrapperVectorFunction(2, 2, f_gradient);
  size_t dim = 2;

  //Create the grid object.
  std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createBsplineGrid(dim, 2));
  grid->getGenerator().regular(4);

  //Sample the gradient function at the grid points
  sgpp::base::GridSample<sgpp::base::DataVector> sample =
      sgpp::base::SampleHelper::sampleGrid(grid, funcGradient);

  //Create the reducer
  auto reducer = sgpp::base::AsQuadReducer();


  //// Use the reducer to first evaluate the function using the sampled gradients
  //sgpp::base::AsInfo info = reducer.evaluate(sample);

  //// Print out all the information that the reducer has gathered
  //for (size_t d = 0; d < dim; ++d) {
  //  std::cout << "dimension: " + std::to_string(d) << std::endl;
  //
  //  sgpp::base::DataVector e(dim);
  //  info.eigenVectors.getColumn(d, e);
  //  std::cout << "eigen vector: " + e.toString() << std::endl;

  //  std::cout << "eigen value: " + std::to_string(info.eigenValues[d]) << std::endl;
  //  std::cout << std::endl;
  //}


  ////Create the cutter used to cut off some dimensions
  ////In this case we cut every dimension that has an eigen value of less than 0.1
  //auto cutter = sgpp::base::AsQuadEigenValueCutter(0.1);
  //// Alternatively, we can also reduce the dimensions to a fixed parameter without caring about the accuracy of that reduction.
  //// auto cutter = sgpp::base::AsQuadFixedCutter(1);

  //
  //// Use the cutter to remove unwanted dimensions
  //sgpp::base::AsResult result = cutter.cut(sample, info);
  //std::cout << "transformation matrix used to reduce the function: " +
  //                 result.transformation.toString()
  //          << std::endl;

  // Apply the result to the function
  // In this step, the function input is transformed using the transformation matrix from the result.
  //auto reducedFunc = result.apply(func);

  //std::cout << reducedFunc->eval(sgpp::base::DataVector(1, 0.4)) << std::endl;
  return 0;
}
