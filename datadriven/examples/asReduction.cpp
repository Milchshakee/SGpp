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
#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/SumFunction.hpp>

size_t dims = 8;
size_t samples = 1000;
size_t maxLevel = 3;

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

double ebolaFunc(const sgpp::base::DataVector& x) {
  return (x[0] + ((x[1] * x[3] * x[4]) / x[6]) + (x[2] * x[7] / x[5])) / (x[4] + x[7]);
}

double borehole(const sgpp::base::DataVector& x) {
  return (2 * M_PI * x[2] * (x[3] - x[5])) /
         (std::log(x[1] / x[0]) *
          (1 + (2 * x[6] * x[2] / (std::log(x[1] / x[0]) * x[0] * x[0] * x[7])) + (x[2] / x[4])));
}



int main(int argc, char* argv[]) {
  sgpp::base::DataMatrix m(2, maxLevel);

  std::shared_ptr<sgpp::base::ScalarFunction> func =
      std::make_shared<sgpp::base::WrapperScalarFunction>(dims, borehole);
  sgpp::base::BoundingBox liberiaBb({{0.1, 0.4},        // beta_1
                                     {0.1, 0.4},        // beta_2
                                     {0.05, 0.2},       // beta_3
                                     {0.41, 1},         // p_1
                                     {0.0276, 0.1702},  // gamma_1
                                     {0.081, 0.21},     // gamma_2
                                     {0.25, 0.5},       // omega
                                     {0.0833, 0.7}});   // psi

  
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
  auto distSample = sgpp::base::DistributionSample(samples, dist);

  for (size_t l = 0; l < maxLevel; l++) {
    std::shared_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dims));
    grid->getGenerator().regular(l);
    sgpp::base::SGridSample sample(grid, *unitFunc);
    sample.hierarchise();
    auto eval = std::make_shared<sgpp::base::InterpolantScalarFunction>(sample);
    auto funcs = std::vector<std::shared_ptr<sgpp::base::ScalarFunction>>{unitFunc, eval};
    auto errorFunc = std::make_shared<sgpp::base::SumFunction>(funcs, std::vector<bool>{true, false});
    double error = sgpp::base::DimReduction::calculateMcL2Error(*errorFunc, distSample);
    double totalError = sgpp::base::DimReduction::calculateMcL2Error(*unitFunc, distSample);
    double relError = (error / totalError);

    m.set(0, l, grid->getSize());
    m.set(1, l, relError);
    }
  m.transpose();
  m.toFile("borehole-org.txt");
}
