#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "sgpp/base/tools/dimension/AnovaReducer.hpp"
#include "sgpp/optimization/function/scalar/WrapperScalarFunction.hpp"

double f(const sgpp::base::DataVector& v) { return v[0]; }

int main(int argc, char* argv[]) {
  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
  size_t dim = 2;
  std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createAnovaBoundaryGrid(dim));
  grid->getGenerator().regular(4);

  sgpp::base::GridSample<double> sample = sgpp::base::SampleHelper::sampleGrid(grid, func);

  auto reducer = sgpp::base::AnovaReducer();
  sgpp::base::AnovaInfo info = reducer.evaluate(sample);
  sgpp::base::AnovaVarianceCutter cutter(0.0);

  sgpp::base::AnovaResult result = cutter.cut(sample, info);
  auto reducedSample = result.apply(sample);

  std::cout << reducedSample.getDimensions() << std::endl;
}
