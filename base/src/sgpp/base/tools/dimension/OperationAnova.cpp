#include <sgpp/base/tools/dimension//OperationAnova.hpp>

sgpp::base::Sample<sgpp::base::AnovaBoundaryGrid::AnovaComponent, double>
sgpp::base::OperationAnova::calculateAnovaComponentVariances(const DataVector& alpha, ErrorRule& rule) {
  std::map<AnovaBoundaryGrid::AnovaComponent, double> variances;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    GridPoint& gp = gridStorage.getPoint(i);
    AnovaBoundaryGrid::AnovaComponent c = AnovaBoundaryGrid::getAnovaComponentOfPoint(gp);
    double val = rule.calculateBasisFunction(alpha[i], gp);
    variances.emplace(c, variances.find(c) == variances.end() ? val : variances.at(c) + val);
  }
  return Sample<AnovaBoundaryGrid::AnovaComponent, double>(variances);
}
