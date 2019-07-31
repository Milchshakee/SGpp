#include <sgpp/base/tools/dimension//OperationAnova.hpp>

double getL2NormOfBasis(const sgpp::base::GridPoint& gp) {
  size_t levelSum = 0;
  size_t level1Count = 0;
  for (size_t d = 0; d < gp.getDimension(); d++) {
    levelSum += std::max<size_t>(gp.getLevel(d), 1) - 1;
    if (gp.getLevel(d) == 1) {
      level1Count++;
    }
  }

  double f = std::pow(2.0 / 3.0, static_cast<double>(gp.getDimension()));
  double exp = std::pow(2.0, -(static_cast<double>(levelSum) - static_cast<double>(level1Count)));
  double val = std::sqrt(f * exp);

  return val;
}

sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double>
sgpp::base::OperationAnova::calculateAnovaComponentVariances(const DataVector& alpha) {
  std::map<AnovaBoundaryGrid::AnovaComponent, double> variances;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    GridPoint& gp = gridStorage.getPoint(i);
    AnovaBoundaryGrid::AnovaComponent c = AnovaBoundaryGrid::getAnovaComponentOfPoint(gp);
    double integral = getL2NormOfBasis(gp);
    double val = std::abs(alpha[i]) * integral;
    variances.emplace(c, variances.find(c) == variances.end() ? val : variances.at(c) + val);
  }
  return Sample<AnovaBoundaryGrid::AnovaComponent, double>(variances);
}
