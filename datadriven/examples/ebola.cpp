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

size_t dims = 8;
size_t samples = 1000;
size_t reducedDims = 5;

double ebolaFunc(const sgpp::base::DataVector& x) {
  return (x[0] + ((x[1] * x[3] * x[4]) / x[6]) + (x[2] * x[7] / x[5])) / (x[4] + x[7]);
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
  std::vector<sgpp::base::DistributionType> types(dims, sgpp::base::DistributionType::Uniform);
  std::shared_ptr<sgpp::base::VectorFunction> bbTrans =
      std::make_shared<sgpp::base::BoundingBoxFunction>(
          sgpp::base::BoundingBoxFunction::Type::FROM_UNIT_BB, liberiaBb);
  auto v = {bbTrans};
  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc =
      std::make_shared<sgpp::base::ChainScalarFunction>(v, func);

  sgpp::datadriven::ASMatrixGradientMC as(unitFunc);
  as.createMatrixMonteCarloFiniteDifference(samples);
  as.evDecompositionForSymmetricMatrices();
  auto a = as.getEigenvaluesDataVector();
    sgpp::base::DistributionsVector dist(types, liberiaBb, true);
  auto distSample = sgpp::base::DistributionSample(samples, dist);
  sgpp::base::PointSample<double> sample = sgpp::base::DimReduction::createActiveSubspaceSample(
      sgpp::base::SampleHelper::sampleScalarFunction(distSample, *unitFunc), as.getTransformationMatrixDataMatrix(reducedDims), reducedDims);
  sgpp::datadriven::Dataset data = sgpp::base::SampleHelper::fromPointSample(sample);

  sgpp::base::DimReduction::RegressionConfig config (3);
  config.gridLevel = 5;
  config.maxIterations = 100;
  config.samples = 100;
  config.regularizationBases = {0.5};
  config.lambdas = {0.1, 0.5, 1};
  sgpp::base::AsReductionResult result = sgpp::base::DimReduction::reduceAS(
      unitFunc, dist, {config, config, config, config, config, config, config, config, config, config});
  double error = sgpp::base::DimReduction::calculateMcL2Error(*result.errorFunction, distSample);
  double totalVar = sgpp::base::DimReduction::calculateMcL2Error(*unitFunc, distSample);
  double va = 0;

  sgpp::base::DistributionsVector dist2(types, liberiaBb, true);
  auto distSample2 = sgpp::base::DistributionSample(10000, dist);
  sgpp::base::PointSample<double> sample2 = sgpp::base::DimReduction::createActiveSubspaceSample(
      sgpp::base::SampleHelper::sampleScalarFunction(distSample2, *unitFunc),
      as.getTransformationMatrixDataMatrix(reducedDims), reducedDims);
  sgpp::datadriven::Dataset data2 = sgpp::base::SampleHelper::fromPointSample(sample2);
  double aa = 0;
}
