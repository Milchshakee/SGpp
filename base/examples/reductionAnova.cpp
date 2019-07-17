#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "sgpp/base/tools/dimension/AnovaReducer.hpp"
#include "sgpp/optimization/function/scalar/WrapperScalarFunction.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/operation/hash/OperationEvalAnovaBoundary.hpp"
#include <iostream>

double f(const sgpp::base::DataVector& v) {  return 2.0 + v[0] + v[1]; }

int main(int argc, char* argv[]) {
  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
  size_t dim = 2;
  std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createAnovaBoundaryGrid(dim));
  grid->getGenerator().regular(2);

  sgpp::base::GridSample<double> sample = sgpp::base::SampleHelper::sampleGrid(grid, func);
  sample.hierarchise();

  //std::unique_ptr<sgpp::base::OperationEval> op(
  //    sgpp::op_factory::createOperationEval(*grid));
  //sgpp::base::DataVector p = sgpp::base::DataVector(2);
  //p[0] = 0.5;
  //p[1] = 0.5;
  //std::cout << op->eval(alpha, p) << std::endl;

  auto reducer = sgpp::base::AnovaReducer();
  sgpp::base::AnovaInfo info = reducer.evaluate(sample);

  sgpp::base::AnovaDimensionVarianceShareCutter cutter(0.95);
  // sgpp::base::AnovaComponentVarianceCutter cutter(0.5);

  
  //sgpp::base::AnovaFixedCutter cutter(1);

  sgpp::base::AnovaResult result = cutter.cut(sample, info);

  
    std::cout << "covered variance: " + std::to_string(result.coveredVariance) << std::endl;

  // Loop over all dimensions and check if they are removed from the reduced grid
  for (size_t d = 0; d < result.dimensions; ++d) {
    std::cout << "dimension: " + std::to_string(d) << std::endl;
    std::cout << "active: " + std::to_string(result.activeDimensions[d]) << std::endl;
    }

  // Loop over all active ANOVA components of the grid
  for (sgpp::base::AnovaHelper::AnovaComponent component : result.activeComponents) {
      std::string str(component.begin(), component.end());
      std::cout << "component: " + str << std::endl;
      std::cout << "component variance: " + std::to_string(info.variances.getValue(component)) << std::endl;
  }


  auto reducedSample = result.apply(sample);

  std::cout << reducedSample.getDimensions() << std::endl;
}
