#include <iostream>

#include "sgpp/base/datatypes/DataMatrix.hpp"
#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include "sgpp/base/tools/dimension/AsMcReducer.hpp"
#include "sgpp/optimization/function/vector/WrapperVectorFunction.hpp"
#include "sgpp/optimization/function/scalar/WrapperScalarFunction.hpp"
#include "sgpp/base/operation/hash/OperationQuadrature.hpp"
#include "sgpp/base/tools/dimension/AsQuadReducer.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"

double f(const sgpp::base::DataVector& v)
{ return v[0]; }

void f_gradient(const sgpp::base::DataVector& x, sgpp::base::DataVector& out) {
  out[0] = 1;
  out[1] = 0;
}

int main() {
//  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
//  auto dist = sgpp::base::RandomUniformDistribution(100, std::mt19937_64::default_seed, 2);
//  auto funcGradient = sgpp::optimization::WrapperVectorFunction(2, 2, f_gradient);
//  sgpp::base::Sample<sgpp::base::DataVector> gradientSamples =
//      sgpp::base::Sampler::sampleVectorFunction(dist, funcGradient);
//  sgpp::base::Sample<sgpp::base::DataMatrix> m =
//      sgpp::base::AsMcReducer::fromGradientSample(gradientSamples);
//
//  auto cutoff = sgpp::base::FixedCutter<sgpp::base::McActiveSubspaceInfo>(1);
//  auto reducer = sgpp::base::AsMcReducer();
//  auto i = sgpp::base::McActiveSubspaceInfo();
//
//  auto result = reducer.evaluateAndReduce(m, cutoff, i);
//  auto reducedFunc = result->apply(func);
//
//  std::cout << reducedFunc->eval(sgpp::base::DataVector(1, 0.4)) << std::endl;
//}
//
//int main2() {
//  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
//  auto dist = sgpp::base::RandomUniformDistribution(100, std::mt19937_64::default_seed, 2);
//  auto funcGradient = sgpp::optimization::WrapperVectorFunction(2, 2, f_gradient);
//  sgpp::base::Sample<sgpp::base::DataMatrix> m =
//      sgpp::base::AsMcReducer::fromFiniteDifferences(func, dist);
//
//  auto cutoff = sgpp::base::FixedCutter<sgpp::base::McActiveSubspaceInfo>(1);
//  auto reducer = sgpp::base::AsMcReducer();
//  auto i = sgpp::base::McActiveSubspaceInfo();
//
//  auto result = reducer.evaluateAndReduce(m, cutoff, i);
//  auto reducedFunc = result->apply(func);
//
//  std::cout << reducedFunc->eval(sgpp::base::DataVector(1, 0.4)) << std::endl;
//}
//
//int main3() {
//  auto func = sgpp::optimization::WrapperScalarFunction(2, f);
//  auto funcGradient = sgpp::optimization::WrapperVectorFunction(2, 2, f_gradient);
//  size_t dim = 2;
//  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
//  grid->getGenerator().regular(4);
//
//  sgpp::base::Sample<sgpp::base::DataVector> sample =
//      sgpp::base::Sampler::sampleGrid(*grid, funcGradient);
//  sgpp::base::Sample<sgpp::base::DataMatrix> m =
//      sgpp::base::AsMcReducer::fromGradientSample(sample);
//
//  
//  auto cutoff = sgpp::base::FixedCutter<sgpp::base::McActiveSubspaceInfo>(1);
//  auto reducer = sgpp::base::AsMcReducer();
//  auto i = sgpp::base::McActiveSubspaceInfo();
//
//  auto result = reducer.evaluateAndReduce(m, cutoff, i);
//  auto reducedFunc = result->apply(func);
//
//  std::cout << reducedFunc->eval(sgpp::base::DataVector(1, 0.4)) << std::endl;
  return 0;
}
