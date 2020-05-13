#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>
#include <sgpp/base/tools/DistributionSample.hpp>
#include <sgpp/base/tools/Sample.hpp>
#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/base/function/vector/BoundingBoxFunction.hpp>
#include <sgpp/base/function/scalar/ChainScalarFunction.hpp>

size_t dims = 8;
size_t samples = 10000;

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
  sgpp::base::DistributionsVector dist(types, liberiaBb);
  std::shared_ptr<sgpp::base::VectorFunction> bbTrans =
      std::make_shared<sgpp::base::BoundingBoxFunction>(
          sgpp::base::BoundingBoxFunction::Type::FROM_UNIT_BB, liberiaBb);
  auto v = {bbTrans};
  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc =
      std::make_shared<sgpp::base::ChainScalarFunction>(v, func);

  sgpp::datadriven::ASMatrixGradientMC as(unitFunc);
  as.createMatrixMonteCarloFiniteDifference(samples);
  as.evDecompositionForSymmetricMatrices();

  auto distSample = sgpp::base::DistributionSample(10000, dist);
  sgpp::base::PointSample<double> sample = sgpp::base::DimReduction::createActiveSubspaceSample(
      sgpp::base::SampleHelper::sampleScalarFunction(distSample, *unitFunc), as.getEigenvectorsDataMatrix(), 1);
  sgpp::datadriven::Dataset data = sgpp::base::SampleHelper::fromPointSample(sample);
}
