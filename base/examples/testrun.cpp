#include <iostream>
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/operation/hash/OperationEvalAnovaBoundary.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include "sgpp/base/tools/dimension/AnovaReducer.hpp"
#include "sgpp/base/function/scalar/WrapperScalarFunction.hpp"
#include "sgpp/base/tools/dimension/AsMcReducer.hpp"
#include "sgpp/base/function/vector/WrapperVectorFunction.hpp"
#include "sgpp/base/operation/hash/OperationHierarchisationAnovaBoundary.hpp"
#include "sgpp/base/tools/dimension/AsQuadReducer.hpp"

double f(const sgpp::base::DataVector& v) {
  return 2.0 + v[0] + v[1];
  //return 2.0 + std::pow(v[0], 5) * std::exp(-v[1]);
}

void f_gradient(const sgpp::base::DataVector& x, sgpp::base::DataVector& out) {
  out[0] = 1;
  out[1] = 1;
  //out[0] = 5 * std::pow(x[0], 4) * std::exp(-x[1]);
  //out[1] = std::pow(x[0], 5) * (-x[1] * std::exp(-x[1]));
}

size_t dim = 2;
size_t maxLevel = 10;
size_t paths = 10000;

//int main(int argc, char* argv[]) {
//  sgpp::base::OperationHierarchisationAnovaBoundary::setIntegralPolicy();
//  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
//  auto reducer = sgpp::base::AnovaReducer();
//
//  sgpp::base::DataMatrix m1(2, maxLevel);
//  sgpp::base::DataMatrix m2(2, maxLevel);
//  for (size_t l = 0; l < maxLevel; ++l) {
//    std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createAnovaBoundaryGrid(dim));
//    grid->getGenerator().regular(l);
//    m1.set(0, l, grid->getStorage().getSize());
//    sgpp::base::SGridSample sample(grid, func);
//    sample.hierarchise();
//    double normalL2Error = sample.mcL2Error(func, paths);
//    m1.set(1, l, normalL2Error);
//
//    sgpp::base::AnovaInfo info = reducer.evaluate(sample);
//    sgpp::base::AnovaFixedCutter cutter(1);
//    sgpp::base::AnovaResult result = cutter.cut(sample, info);
//    auto& reducedSample = result.getReducedOutput();
//    m2.set(0, l, reducedSample.getGrid().getSize());
//    //double test = reducedSample.eval(sgpp::base::DataVector(1, 0.4));
//    double reducedL2Error = result.calcMcL2Error(func, paths);
//    m2.set(1, l, reducedL2Error);
//  }
//  m1.transpose();
//  m1.toFile("results1.txt");
//  m2.transpose();
//  m2.toFile("results2.txt");
//}

int main(int argc, char* argv[]) {
  auto func = sgpp::base::WrapperScalarFunction(2, f);
  auto reducer = sgpp::base::AsQuadReducer();

  sgpp::base::DataMatrix m1(2, maxLevel);
  sgpp::base::DataMatrix m2(2, maxLevel);
  for (size_t l = 0; l < maxLevel; ++l) {
    std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
    grid->getGenerator().regular(l);
    m1.set(0, l, grid->getStorage().getSize());
    sgpp::base::SGridSample sample(grid, func);
    sample.hierarchise();
    double normalL2Error = sample.mcL2Error(func, paths);
    m1.set(1, l, normalL2Error);
    sgpp::base::EvalFunction fe(sample);

    auto dist = sgpp::base::RandomUniformDistribution(paths, std::mt19937_64::default_seed, 2);
      auto funcGradient = sgpp::base::WrapperVectorFunction(2, 2, f_gradient);
      sgpp::base::GridSample<sgpp::base::DataVector> gradientSamples = sgpp::base::SampleHelper::sampleGrid(grid, funcGradient);

    sgpp::base::AsQuadInput i{sample, gradientSamples};
    sgpp::base::AsInfo info = reducer.evaluate(i);
      sgpp::base::AsQuadFixedCutter cutter(1);
    sgpp::base::AsQuadResult result = cutter.cut(i, info);
    auto& reducedSample = result.getReducedOutput();
    // double test = reducedSample.eval(sgpp::base::DataVector(1, 0.4));
    double reducedL2Error = result.calcMcL2Error(func, paths);

    
    std::shared_ptr<sgpp::base::Grid> grid2(sgpp::base::Grid::createLinearBoundaryGrid(dim - 1));
    grid2->getGenerator().regular(l);
    m2.set(0, l, grid2->getStorage().getSize());

    m2.set(1, l, reducedL2Error);
  }
  m1.transpose();
  m2.transpose();
  m1.toFile("results1.txt");
  m2.toFile("results2.txt");
}

//int main(int argc, char* argv[]) {
//  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
//  auto reducer = sgpp::base::AsMcReducer();
//
//  sgpp::base::DataMatrix m1(2, maxLevel);
//  sgpp::base::DataMatrix m2(2, maxLevel);
//  for (size_t l = 0; l < maxLevel; ++l) {
//    std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
//    grid->getGenerator().regular(l);
//    m1.set(0, l, grid->getStorage().getSize());
//    sgpp::base::SGridSample sample(grid, func);
//    sample.hierarchise();
//    double normalL2Error = sample.mcL2Error(func, paths);
//    m1.set(1, l, normalL2Error);
//    sgpp::base::EvalFunction fe(sample);
//
//    auto dist = sgpp::base::RandomUniformDistribution(paths, std::mt19937_64::default_seed, 2);
//    auto funcGradient = sgpp::optimization::WrapperVectorFunction(2, 2, f_gradient);
//    sgpp::base::PointSample<sgpp::base::DataVector> gradientSamples =
//        sgpp::base::SampleHelper::sampleVectorFunction(dist, funcGradient);
//    sgpp::base::PointSample<sgpp::base::DataMatrix> m =
//        sgpp::base::AsMcReducer::fromGradientSample(gradientSamples);
//    sgpp::base::AsMcInput i{fe, m};
//
//    sgpp::base::AsInfo info = reducer.evaluate(i);
//    sgpp::base::AsMcFixedCutter cutter(1);
//    sgpp::base::AsMcResult result = cutter.cut(i, info);
//    auto& reducedSample = result.getReducedOutput();
//    // double test = reducedSample.eval(sgpp::base::DataVector(1, 0.4));
//    double reducedL2Error = result.calcMcL2Error(func, paths);
//
//    std::shared_ptr<sgpp::base::Grid> grid2(sgpp::base::Grid::createLinearBoundaryGrid(dim - 1));
//    grid2->getGenerator().regular(l);
//    m2.set(0, l, grid2->getStorage().getSize());
//
//    m2.set(1, l, reducedL2Error);
//  }
//  m1.transpose();
//  m2.transpose();
//  m1.toFile("results1.txt");
//  m2.toFile("results2.txt");
//}
