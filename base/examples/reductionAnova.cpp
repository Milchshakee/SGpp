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
  // Create anchor at (0,0)
  std::vector<sgpp::base::AnovaTypes::LevelIndexPair> anchor(2, {-1, 0});
  std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createAnovaLinearBoundaryGrid(dim));
  grid->getGenerator().regular(2);

  sgpp::base::SGridSample sample(grid, func);
  sample.hierarchise();
  sgpp::base::AnovaInput input{sample, func};

  sgpp::base::L2SquaredMcRule rule(std::mt19937_64::default_seed, 10000);

  auto reducer = sgpp::base::AnovaReducer(rule);
  sgpp::base::AnovaInfo info = reducer.evaluate(input);

  // Preserve a share of 0.95 of the total variance
  //sgpp::base::AnovaErrorRuleCutter cutter(rule, 0.9);

  sgpp::base::AnovaFixedCutter cutter(1);

  sgpp::base::AnovaResult result = cutter.cut(input, info);

  
    std::cout << "covered variance: " + std::to_string(result.coveredVariance) << std::endl;

  // Loop over all dimensions and check if they are removed from the reduced grid
  for (size_t d = 0; d < sample.getDimensions(); ++d) {
    std::cout << "dimension: " + std::to_string(d) << std::endl;
    std::cout << "active: " + std::to_string(result.activeDimensions[d]) << std::endl;
    }

  // Loop over all ANOVA components of the grid
  for (size_t i = 0; i < info.variances.getSize(); ++i) {
      sgpp::base::AnovaBoundaryGrid::AnovaComponent& component = info.variances.getKeys()[i];
    std::string str(component.begin(), component.end());
    std::cout << "component: " + str << std::endl;
    std::cout << "variance component share: " + std::to_string(info.variances.getValues()[i])
              << std::endl;
  }

  sgpp::base::SGridSample& reducedSample = result.getReducedOutput();

  std::cout << "reduced dimensions: " << reducedSample.getDimensions() << std::endl;

  sgpp::base::L2McRule l2Rule(std::mt19937_64::default_seed, 10000);
  double l2Error = std::sqrt(result.calculateAbsoluteError(l2Rule));
  std::cout << "L2 error: " << l2Error << std::endl;
}
