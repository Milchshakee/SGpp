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

sgpp::base::AsReductionResult reduceBorehole() {
  size_t dims = 8;
  std::shared_ptr<sgpp::base::ScalarFunction> func =
      std::make_shared<sgpp::base::WrapperScalarFunction>(dims, borehole);
  sgpp::base::BoundingBox boreholeBb({{0.05, 0.15},     // r_w
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
          sgpp::base::BoundingBoxFunction::Type::FROM_UNIT_BB, boreholeBb);
  auto v = {bbTrans};
  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc =
      std::make_shared<sgpp::base::ChainScalarFunction>(v, func);

  sgpp::base::DistributionsVector dist(types, boreholeBb, true, chars);

  sgpp::base::DimReduction::RegressionConfig config(2);
  config.gridLevel = 3;
  config.maxIterations = 5000;
  config.samples = 1000;
  // config.regularizationBases = {1.0, 0.5};
  config.lambdas = {0.1, 0.5, 1};
  config.trainDataShare = 0.02;
  config.refinements = 1;
  config.refinementPoints = 100;
  config.crossValidations = 10;
  config.subIterations = 10;

  sgpp::base::DimReduction::RegressionConfig c1 = config;
  c1.reducedDimension = 4;
  c1.gridLevel = 1;
  c1.subIterations = 5;

  sgpp::base::DimReduction::RegressionConfig c2 = config;
  c2.reducedDimension = 6;
  c2.gridLevel = 0;
  c2.subIterations = 3;

  auto configs = std::vector<sgpp::base::DimReduction::RegressionConfig>(iterations, config);
  sgpp::base::DataMatrix mat(2, iterations);
  // configs[0] = c1;
  // configs[1] = c1;
  // configs[2] = c1;
  //configs[6] = c1;
  //configs[7] = c1;
  //configs[8] = c1;
  //configs[9] = c1;
  sgpp::base::AsReductionResult result =
      sgpp::base::DimReduction::reduceAS(unitFunc, dist, samples, configs);
  return result;
}

  double analytical(const sgpp::base::DataVector& x) {
  return std::exp(x[0] + x[1] + x[2] + x[3] + x[4])
         + 5 * (std::sin(M_PI * x[5] * x[6] * x[7]) * x[8] * x[8]) * std::pow(x[9], 5) * x[10] +
           std::pow(x[11] + 2, 2) * x[12] * x[13] * x[14] + x[15] * x[16] +
           std::log(1 + (10 * x[17] / (x[18] + x[19])));
  }

  double gfunction(const sgpp::base::DataVector& x) {
    size_t d = 2;
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

double friedman1(const sgpp::base::DataVector& x) {
    std::normal_distribution<double> dist;
  std::default_random_engine rand;
  return 10.0 * sin(M_PI * x[0] * x[1]) + 20.0 * (x[2] - 0.5) * (x[2] - 0.5) + 10.0 * x[3] +
           5.0 * x[4] + dist(rand);
  }

  sgpp::base::AsReductionResult reduceAnalytical()
{
  size_t dims = 10;
    std::shared_ptr<sgpp::base::ScalarFunction> func =
        std::make_shared<sgpp::base::WrapperScalarFunction>(dims, friedman1);

  std::vector<sgpp::base::DistributionType> types(dims, sgpp::base::DistributionType::Uniform);
  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc = func;

  sgpp::base::BoundingBox bb(dims);
  sgpp::base::DistributionsVector dist(types, bb, true);

  sgpp::base::DimReduction::RegressionConfig config(3);
  config.gridLevel = 1;
  config.maxIterations = 1000;
  config.samples = 1000;
  // config.regularizationBases = {1.0, 0.5};
  //config.lambdas = {0.1, 0.5, 1};
  config.trainDataShare = 0.02;
  config.refinements = 0;
  config.refinementPoints = 5;
  config.crossValidations = 10;
  config.subIterations = 20;

  sgpp::base::DimReduction::RegressionConfig c1 = config;
  c1.reducedDimension = 4;
  c1.gridLevel = 1;
  c1.subIterations = 3;

  sgpp::base::DimReduction::RegressionConfig c2 = config;
  c2.reducedDimension = 6;
  c2.gridLevel = 0;
  c2.subIterations = 3;

  auto configs = std::vector<sgpp::base::DimReduction::RegressionConfig>(iterations, config);
  sgpp::base::DataMatrix mat(2, iterations);
  // configs[2] = c1;
  sgpp::base::AsReductionResult result =
      sgpp::base::DimReduction::reduceAS(unitFunc, dist, samples, configs, false);
  return result;
}

int main()
{
  sgpp::base::DataMatrix mat(2, iterations);
  sgpp::base::AsReductionResult result = reduceAnalytical();

    for (size_t i = 0; i < iterations; i++) {
    size_t sum = 0;
      for (size_t j = 0; j <= i; j++) {
      sum += result.reductions[j].gridPoints;
      }
      mat(0, i) = sum;
      mat(1, i) = result.reductions[i].l2Error;
  }

  mat.transpose();
  mat.toFile("friedman1.txt");
}
