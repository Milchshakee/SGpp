#include <sgpp/base/tools/dimension//OperationAnova.hpp>

sgpp::base::Sample<sgpp::base::AnovaHelper::AnovaComponent, double> sgpp::base::OperationAnova::calculateAnovaComponentVariances(
    const DataVector& alpha) {
  return AnovaHelper::calculateAnovaComponentVariances(gridStorage, alpha);
}
