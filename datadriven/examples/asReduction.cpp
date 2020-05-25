#include <sgpp/base/function/scalar/ChainScalarFunction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/function/vector/BoundingBoxFunction.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/DistributionSample.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>
#include <sgpp/base/tools/EigenHelper.hpp>
#include <sgpp/base/tools/Sample.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/datadriven/tools/dimension/DimReduction.hpp>

size_t dims = 2;
size_t samples = 1000;
size_t reducedDims = 1;

double sinfunc(const sgpp::base::DataVector& x) {
  return x[1];
}

double fsin = std::sin(0.25 * M_PI);
double fcos = std::cos(0.25 * M_PI);

double g(const sgpp::base::DataVector& v) {
  for (size_t t = 0; t < dims; t++) {
    if ((v[t] < 0.0) || (v[t] > 1.0)) {
      return std::numeric_limits<double>::infinity();
    }
  }

  double x = fcos * v[0] - fsin * v[1];
  double y = fsin * v[0] + fcos * v[1];
  return std::max(1 - std::abs(5 * x - 2.5), 0.0);
}

int main(int argc, char* argv[]) {
  std::shared_ptr<sgpp::base::ScalarFunction> func =
      std::make_shared<sgpp::base::WrapperScalarFunction>(dims, g);
  std::vector<sgpp::base::DistributionType> types(dims, sgpp::base::DistributionType::Uniform);

  auto bb = sgpp::base::BoundingBox(dims);
  sgpp::base::DistributionsVector dist(types, bb, true);

  sgpp::base::DimReduction::RegressionConfig config(1);
  config.gridLevel = 10;
  config.samples = 1000;
  auto c2 = config;
  c2.reducedDimension = 2;
  sgpp::base::AsReductionResult result = sgpp::base::DimReduction::reduceAS(
      func, dist,
      {config, c2});

  auto distSample = sgpp::base::DistributionSample(samples, dist);
    double error = sgpp::base::DimReduction::calculateMcL2Error(*result.errorFunction, distSample);
  double totalVar = sgpp::base::DimReduction::calculateMcL2Error(*func, distSample);
  sgpp::base::DataVector x;

  sgpp::base::DataVector v = sgpp::base::DataVector(2, 1);
  v[1] = 0;
   result.reductions[0].transformation->eval(v, x);
  double va = 0;
}
