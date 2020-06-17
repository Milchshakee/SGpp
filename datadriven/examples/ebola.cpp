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

size_t iterations = 5;
size_t samples = 10000;

double ebolaFunc(const sgpp::base::DataVector& x) {
  return (x[0] + ((x[1] * x[3] * x[4]) / x[6]) + (x[2] * x[7] / x[5])) / (x[4] + x[7]);
}

double borehole(const sgpp::base::DataVector& x) {
  return (2 * M_PI * x[2] * (x[3] - x[5])) / (std::log(x[1] / x[0]) * (1 + (2 * x[6] * x[2] / (std::log(x[1] / x[0]) * x[0] * x[0] * x[7])) + (x[2] / x[4])));
}

double test1(const sgpp::base::DataVector& x) { return std::sin(M_PI * x[0]) + x[1]; }

double friedman1(const sgpp::base::DataVector& x) {
  std::normal_distribution<double> dist;
  std::default_random_engine rand;
  return 10.0 * sin(M_PI * x[0] * x[1]) + 20.0 * (x[2] - 0.5) * (x[2] - 0.5) + 10.0 * x[3] +
         5.0 * x[4] + dist(rand);
}

//sgpp::base::AsReductionResult reduceEbola() {
//  size_t dims = 8;
//  std::shared_ptr<sgpp::base::ScalarFunction> func =
//      std::make_shared<sgpp::base::WrapperScalarFunction>(dims, ebolaFunc);
//  std::shared_ptr<sgpp::base::BoundingBox> liberiaBb = std::make_shared<sgpp::base::BoundingBox>(std::vector<sgpp::base::BoundingBox1D>{
//                                     {0.1, 0.4},        // beta_1
//                                     {0.1, 0.4},        // beta_2
//                                     {0.05, 0.2},       // beta_3
//                                     {0.41, 1},         // p_1
//                                     {0.0276, 0.1702},  // gamma_1
//                                     {0.081, 0.21},     // gamma_2
//                                     {0.25, 0.5},       // omega
//                                     {0.0833, 0.7}});   // psi
//
//  std::vector<sgpp::base::DistributionType> types(dims, sgpp::base::DistributionType::Uniform);
//  std::shared_ptr<sgpp::base::VectorFunction> bbTrans =
//      std::make_shared<sgpp::base::BoundingBoxFunction>(
//          sgpp::base::BoundingBoxFunction::Type::FROM_UNIT_BB, *liberiaBb);
//  auto v = {bbTrans};
//  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc =
//      std::make_shared<sgpp::base::ChainScalarFunction>(v, func);
//
//  sgpp::base::DistributionsVector dist(types, {} ,liberiaBb);
//
//  sgpp::base::DimReduction::RegressionConfig config(1);
//  config.gridLevel = 2;
//  config.maxIterations = 1000;
//  config.samples = 1000;
//  // config.regularizationBases = {1.0, 0.5};
//  config.lambdas = {0.1, 0.5, 1};
//  config.trainDataShare = 0.02;
//  config.refinements = 0;
//  config.refinementPoints = 100;
//  config.crossValidations = 5;
//  config.subIterations = 10;
//
//  sgpp::base::DimReduction::RegressionConfig c1 = config;
//  c1.reducedDimension = 2;
//  c1.gridLevel = 2;
//  c1.subIterations = 5;
//
//  sgpp::base::DimReduction::RegressionConfig c2 = config;
//  c2.reducedDimension = 6;
//  c2.gridLevel = 0;
//  c2.subIterations = 3;
//
//  auto configs = std::vector<sgpp::base::DimReduction::RegressionConfig>(iterations, config);
//  sgpp::base::DataMatrix mat(2, iterations);
//  // configs[7] = c1;
//  // configs[8] = c1;
//  // configs[9] = c1;
//  sgpp::base::DimReduction::reduceASRandom(unitFunc, dist, config, 1);
//  sgpp::base::AsReductionResult result =
//      sgpp::base::DimReduction::reduceAS(unitFunc, dist, samples, configs);
//  return result;
//}

sgpp::base::AsReductionResult reduceBorehole() {
  size_t dims = 2;
  std::shared_ptr<sgpp::base::ScalarFunction> func =
      std::make_shared<sgpp::base::WrapperScalarFunction>(dims, test1);
  std::shared_ptr<sgpp::base::BoundingBox> boreholeBb = std::make_shared<sgpp::base::BoundingBox>(std::vector<sgpp::base::BoundingBox1D>{
                                      {0.05, 0.15},     // r_w
                                      {100, 50000},     // r
                                      {63070, 115600},  // T_u
                                      {990, 1110},      // H_u
                                      {63.1, 116},      // T_I
                                      {700, 820},       // H_I
                                      {1120, 1680},     // L
                                      {9855, 12045}});  // K_W

  std::vector<sgpp::base::DistributionType> types(dims, sgpp::base::DistributionType::Uniform);
  types[0] = sgpp::base::DistributionType::TruncNormal;
  types[1] = sgpp::base::DistributionType::TruncLognormal;
  std::vector<sgpp::base::DataVector> chars(dims);
  chars[0] = sgpp::base::DataVector{0.1, 0.01618};
  chars[1] = sgpp::base::DataVector{7.71, 1.0056};
  std::shared_ptr<sgpp::base::VectorFunction> bbTrans =
      std::make_shared<sgpp::base::BoundingBoxFunction>(
          sgpp::base::BoundingBoxFunction::Type::FROM_UNIT_BB, *boreholeBb);
  auto v = {bbTrans};
  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc =
      std::make_shared<sgpp::base::ChainScalarFunction>(v, func);

  std::shared_ptr<sgpp::base::BoundingBox> bb = std::make_shared<sgpp::base::BoundingBox>(dims);
  sgpp::base::DistributionsVector dist(
      {sgpp::base::DistributionType::Uniform, sgpp::base::DistributionType::Uniform}, {}, bb);
  sgpp::base::DistributionSample distSample(2000, dist);
  sgpp::base::PointSample<double> funcSample =
      sgpp::base::SampleHelper::sampleScalarFunction(distSample, *func);

  sgpp::base::DimReduction::ReductionConfig redConfig(1);

  sgpp::base::DimReduction::RegressionConfig regConfig;
  regConfig.sampleCount = 1000;
  regConfig.maxGridPoints = 100;
  regConfig.refinementPointsPerGridPoint = 0;
  regConfig.trainDataPerGridPoint = 1;
  regConfig.crossValidations = 5;

  sgpp::base::DimReduction::BasisConfig basis;
  basis.type = sgpp::base::DimReduction::BasisConfig::INV_AS;
  basis.basisIterations = 5;
  basis.evaluationConfig = regConfig;

  std::shared_ptr<sgpp::base::VectorFunction> gradient =
      sgpp::base::DimReduction::finiteDifferencesFunction(func, 1e-8);
  sgpp::base::DimReduction::GradientConfig gradConfig;
  gradConfig.sampleCount = 1000;
  gradConfig.maxNorm = 3;
  gradConfig.pNorm = 1;
  gradConfig.type = sgpp::base::DimReduction::GradientConfig::GRADIENT_FUNCTION;
  gradConfig.gradientFunction = gradient;

  auto exConfig = sgpp::base::DimReduction::ExaminationConfig{};
  exConfig.discardWorseIterations = false;

  std::vector<sgpp::base::DimReduction::Output> output;

  sgpp::base::DataMatrix mat(2, iterations);

  sgpp::base::AsReductionResult result =
      sgpp::base::DimReduction::reduceAS(funcSample, redConfig, regConfig, gradConfig, basis, exConfig, output);
  return result;
}

  double analytical(const sgpp::base::DataVector& x) {
  return std::exp(0.2 * (x[0] + x[1] + x[2] + x[3] + x[4]))
         + 5 * (std::sin(M_PI * x[5] * x[6] * x[7]) * x[8] * x[8]) * std::pow(x[9], 5) * x[10] +
           std::pow(x[11] + 2, 2) * x[12] +
      2 * x[13] * 3 * x[14] +
      x[15] * x[16] +
           std::log(1 + (10 * x[17] / (0.1 + x[18] + x[19])));
  }

  double gfunction(const sgpp::base::DataVector& x) {
    size_t d = 10;
    double mult = 1;
    for (size_t j = 0; j < d; j++) {
      double a = ((j + 1) - 2) / 2.0;
      mult *= (std::abs(4 * x[j] - 2) + a) / (1 + a);
    }
    return mult;
  }

double morokoffcaflisch(const sgpp::base::DataVector& x)
{
    size_t d = 10;
    double mult = 1 + std::pow(1.0 / d, d);
    for (size_t j = 0; j < d; j++) {
      mult *= std::pow(x[j], (1 / (j + 1.0)));
    }
  }

  double e(const sgpp::base::DataVector& x) {
    return std::exp(x[0] + x[1] + x[2] + x[3] + x[4]);
  }
//
//  sgpp::base::AsReductionResult reduceAnalytical()
//{
//  size_t dims = 20;
//    std::shared_ptr<sgpp::base::ScalarFunction> func =
//        std::make_shared<sgpp::base::WrapperScalarFunction>(dims, analytical);
//
//  std::vector<sgpp::base::DistributionType> types(dims, sgpp::base::DistributionType::Uniform);
//  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc = func;
//
//  std::shared_ptr<sgpp::base::BoundingBox> bb = std::make_shared<sgpp::base::BoundingBox>(dims);
//  sgpp::base::DistributionsVector dist(types, {}, bb);
//
//  sgpp::base::DimReduction::RegressionConfig config(1);
//  config.gridLevel = 4;
//  config.maxIterations = 1000;
//  config.samples = 1000;
//  // config.regularizationBases = {1.0, 0.5};
//  //config.lambdas = {0.1, 0.5, 1};
//  config.trainDataShare = 0.1;
//  config.refinements = 0;
//  config.refinementPoints = 5;
//  config.crossValidations = 5;
//  config.subIterations = 10;
//
//  auto configs = std::vector<sgpp::base::DimReduction::RegressionConfig>(iterations, config);
//  sgpp::base::DataMatrix mat(2, iterations);
//  // configs[2] = c1;
//  sgpp::base::AsReductionResult result =
//      sgpp::base::DimReduction::reduceAS(unitFunc, dist, samples, configs, true);
//  return result;
//}

int main()
{
  sgpp::base::DataMatrix mat(2, iterations);
  sgpp::base::AsReductionResult result = reduceBorehole();

    for (size_t i = 0; i < iterations; i++) {
    size_t sum = 0;
      for (size_t j = 0; j <= i; j++) {
      sum += 1;  // result.reductions[j].gridPoints;
      }
      mat(0, i) = sum;
      mat(1, i) = 0.5;
      //result.reductions[i].l2Error;
  }

  mat.transpose();
  mat.toFile("bore-d1-l7-0.02.txt");
}
