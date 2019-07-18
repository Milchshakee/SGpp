#include <iostream>
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/operation/hash/OperationEvalAnovaBoundary.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "sgpp/base/tools/dimension/AnovaReducer.hpp"
#include "sgpp/optimization/function/scalar/WrapperScalarFunction.hpp"

double f(const sgpp::base::DataVector& v) { return 2.0 + v[0] + v[1]; }

size_t dim = 2;
size_t maxLevel = 10;
size_t paths = 10000;

int main(int argc, char* argv[]) {
  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
  auto reducer = sgpp::base::AnovaReducer();

  sgpp::base::DataMatrix m(3, maxLevel);
  for (size_t l = 0; l < maxLevel; ++l) {
    std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createAnovaBoundaryGrid(dim));
    grid->getGenerator().regular(l);
    m.set(0, l, grid->getStorage().getSize());
    sgpp::base::SGridSample sample(grid, func);
    sample.hierarchise();
    double normalL2Error = sample.mcL2Error(func, paths);
    m.set(1, l, normalL2Error);

    sgpp::base::AnovaInfo info = reducer.evaluate(sample);
    sgpp::base::AnovaDimensionVarianceShareCutter cutter(0.95);
    sgpp::base::AnovaResult result = cutter.cut(sample, info);
    auto& reducedSample = result.getReducedOutput();
    double reducedL2Error = result.calcMcL2Error(func, paths);
    m.set(2, l, normalL2Error);
  }
  m.toFile("results.txt");
}
