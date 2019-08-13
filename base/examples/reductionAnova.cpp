#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/Sample.hpp>
#include <sgpp/base/tools/dimension/AnovaReducer.hpp>
#include <sgpp/base/operation/hash/OperationEvalAnovaLinearBoundary.hpp>
#include <iostream>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>

double f(const sgpp::base::DataVector& v) {  return 2.0 + v[0] + v[1]; }

int main(int argc, char* argv[]) {
  auto func = sgpp::base::WrapperScalarFunction(2, f);
  size_t dim = 2;
  std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createAnovaPrewaveletBoundaryGrid(dim));
  grid->getGenerator().regular(2);

  sgpp::base::SGridSample sample(grid, func);
  sample.hierarchise();

  sgpp::base::L2SquaredMcRule rule(std::mt19937_64::default_seed, 10000);

  auto reducer = sgpp::base::AnovaReducer(rule);
  sgpp::base::AnovaInfo info = reducer.evaluate(sample);

  sgpp::base::AnovaErrorRuleCutter cutter(rule, 0.05);
  //sgpp::base::AnovaFixedCutter cutter(rule, 1);

  sgpp::base::AnovaResult result = cutter.cut(sample, info);

  
    std::cout << "measured relative error: " + std::to_string(result.errorShare) << std::endl;

  // Loop over all dimensions and check if they are removed from the reduced grid
  for (size_t d = 0; d < result.dimensions; ++d) {
    std::cout << "dimension: " + std::to_string(d) << std::endl;
    std::cout << "active: " + std::to_string(result.activeDimensions[d]) << std::endl;
    }

  // Loop over all ANOVA components of the grid
  for (size_t i = 0; i < info.errorShares.getSize(); ++i) {
      sgpp::base::AnovaBoundaryGrid::AnovaComponent& component = info.errorShares.getKeys()[i];
    std::string str(component.begin(), component.end());
    std::cout << "component: " + str << std::endl;
    std::cout << "relative component error: " + std::to_string(info.errorShares.getValues()[i])
              << std::endl;
  }

  sgpp::base::SGridSample& reducedSample = result.getReducedOutput();

  std::cout << "reduced dimensions: " << reducedSample.getDimensions() << std::endl;

  double l2Error = std::sqrt(result.calculateAbsoluteError(rule));
  std::cout << "L2 error: " << l2Error << std::endl;
}
