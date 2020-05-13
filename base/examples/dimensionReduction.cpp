#include <iostream>
#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/globaldef.hpp>
#include "sgpp/base/function/scalar/WrapperScalarFunction.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include <sgpp/base/tools/DistributionUniform.hpp>

double f(const sgpp::base::DataVector& v) { return 3 * v[0] + 1 * v[1]; }

double fsin = std::sin(0.15 * M_PI);
double fcos = std::cos(0.15 * M_PI);

double g(const sgpp::base::DataVector& v) {
  double x = fcos * v[0] - fsin * v[1];
  double y = fsin * v[0] + fcos * v[1];
  return std::max(1 - std::abs(5 * x - 2.5), 0.0);
}


int main(int argc, char* argv[]) {
  sgpp::base::WrapperScalarFunction func(2, g);

  std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(1));
  grid->getGenerator().regular(6);
  std::shared_ptr<sgpp::base::DistributionUniform> u = std::make_shared<sgpp::base::DistributionUniform>();
  sgpp::base::DistributionsVector v(2, u);
  auto dist = sgpp::base::DistributionSample(1000000, v);
  sgpp::base::ActiveSubspaceInfo i = sgpp::base::DimReduction::activeSubspaceMC(func, dist);
  sgpp::base::ReductionResult result = sgpp::base::DimReduction::reduce(func, i.eigenVectors, 1, 7);
  //sgpp::base::EvalFunction eval(result.reducedSample);
  //double va = eval.eval(sgpp::base::DataVector(1, 0.25));
  //sgpp::base::SGridSample b(grid, eval);
  //b.hierarchise();
  //sgpp::base::EvalFunction bEval(b);

  //double err = std::pow(
  //    sgpp::base::DimReduction::calculateMcL2Error(func, result.transformation, bEval, 1, 10000),
  //    2);

  double a = 0;
}