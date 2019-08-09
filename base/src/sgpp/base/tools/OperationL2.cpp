#include <sgpp/base/tools/OperationL2.hpp>
#include <sgpp/base/function/scalar/EvalFunction.hpp>


double sgpp::base::OperationL2::calculateL2Norm(ScalarFunction& func) {
  std::mt19937_64 rand(seed);
  std::uniform_real_distribution<double> dist(0, 1);
  size_t funcDimensions = func.getNumberOfParameters();
  sgpp::base::DataVector point(funcDimensions);
  double res = 0;

  for (size_t i = 0; i < samples; i++) {
    for (size_t d = 0; d < funcDimensions; d++) {
      point[d] = dist(rand);
    }

    double val = func.eval(point);   ;
    res += pow(val, 2);
  }

  return sqrt(res / static_cast<double>(samples));
}


double sgpp::base::OperationL2::calculateMcL2Error(ScalarFunction& func,
  VectorFunction& transformation, ScalarFunction& reduced) {
  std::mt19937_64 rand(seed);
  std::uniform_real_distribution<double> dist(0, 1);
  size_t funcDimensions = func.getNumberOfParameters();
  size_t newDimensions = reduced.getNumberOfParameters();

  sgpp::base::DataVector point(funcDimensions);
  double res = 0;

  for (size_t i = 0; i < samples; i++) {
    for (size_t d = 0; d < funcDimensions; d++) {
      point[d] = dist(rand);
    }
    double val = func.eval(point);
    DataVector out(newDimensions);
    transformation.eval(point, out);
    res += pow(val - reduced.eval(out), 2);
  }

  return sqrt(res / static_cast<double>(samples));
}
