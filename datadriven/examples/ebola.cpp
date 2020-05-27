#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>
#include <sgpp/base/tools/DistributionSample.hpp>
#include <sgpp/base/tools/Sample.hpp>
#include <sgpp/datadriven/tools/dimension/DimReduction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/base/function/vector/BoundingBoxFunction.hpp>
#include <sgpp/base/function/scalar/ChainScalarFunction.hpp>
#include <sgpp/base/tools/EigenHelper.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/SumFunction.hpp>

size_t dims = 8;
size_t samples = 1000;
size_t reducedDims = 5;

double ebolaFunc(const sgpp::base::DataVector& x) {
  return (x[0] + ((x[1] * x[3] * x[4]) / x[6]) + (x[2] * x[7] / x[5])) / (x[4] + x[7]);
}

double borehole(const sgpp::base::DataVector& x) {
  return (2 * M_PI * x[2] * (x[3] - x[5])) / (std::log(x[1] / x[0]) * (1 + (2 * x[6] * x[2] / (std::log(x[1] / x[0]) * x[0] * x[0] * x[7])) + (x[2] / x[4])));
}

int main(int argc, char* argv[]) {
  std::shared_ptr<sgpp::base::ScalarFunction> func = std::make_shared<sgpp::base::WrapperScalarFunction>(dims, ebolaFunc);
  sgpp::base::BoundingBox liberiaBb({
                              {0.1, 0.4}, //beta_1
                              {0.1, 0.4}, //beta_2
                              {0.05, 0.2}, //beta_3
                              {0.41, 1}, //p_1
                              {0.0276, 0.1702}, //gamma_1
                              {0.081, 0.21}, //gamma_2
                              {0.25, 0.5}, //omega
                              {0.0833, 0.7}}); //psi

    sgpp::base::BoundingBox boreholeBb({{0.05, 0.15},        // r_w
                                     {100, 50000},        // r
                                     {63070, 115600},       // T_u
                                     {990, 1110},         // H_u
                                     {63.1, 116},  // T_I
                                     {700, 820},     // H_I
                                     {1120, 1680},       // L
                                     {9855, 12045}});   // K_W

  std::vector<sgpp::base::DistributionType> types(dims, sgpp::base::DistributionType::Uniform);
  std::shared_ptr<sgpp::base::VectorFunction> bbTrans =
      std::make_shared<sgpp::base::BoundingBoxFunction>(
          sgpp::base::BoundingBoxFunction::Type::FROM_UNIT_BB, liberiaBb);
  auto v = {bbTrans};
  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc =
      std::make_shared<sgpp::base::ChainScalarFunction>(v, func);

    sgpp::base::DistributionsVector dist(types, liberiaBb, true);
  auto distSample = sgpp::base::DistributionSample(samples, dist);

  sgpp::base::DimReduction::RegressionConfig config (3);
  config.gridLevel = 2;
  config.maxIterations = 1000;
  config.samples = 1000;
  //config.regularizationBases = {2, 1.25, 1.0, 0.5, 0.25, 0.125};
  config.lambdas = {0.1, 0.5, 1};
  config.trainDataShare = 0.03;
  config.errorCalcSamples = 1000;
  config.refinements = 1;
  config.refinementPoints = 100;

  sgpp::base::DimReduction::RegressionConfig c1(5);
  c1.gridLevel = 1;

  auto configs = std::vector<sgpp::base::DimReduction::RegressionConfig>(20, config);
  sgpp::base::DataMatrix mat(2, 20);
  //configs[0] = c1;
  //configs[1] = c1;
  //configs[2] = c1;
  //configs[3] = c1;
  //configs[4] = c1;
  //configs[5] = c1;
  sgpp::base::AsReductionResult result = sgpp::base::DimReduction::reduceAS(
      unitFunc, dist, configs);

    for (int i = 0; i < 20; i++) {
    size_t sum = 0;
      for (int j = 0; j <= i; j++) {
      sum += result.reductions[j].gridPoints;
      }
      mat(0, i) = sum;;
      mat(1, i) = result.reductions[i].l2Error;
  }

  mat.transpose();
  mat.toFile("ebola-red.txt");

  double error = sgpp::base::DimReduction::calculateMcL2Error(*result.errorFunction, distSample);
  double totalVar = sgpp::base::DimReduction::calculateMcL2Error(*unitFunc, distSample);
  double va = 0;
}
