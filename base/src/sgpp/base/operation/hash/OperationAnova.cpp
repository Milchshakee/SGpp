#include <sgpp/base/operation/hash/OperationAnova.hpp>

sgpp::base::Sample<sgpp::base::AnovaHelper::AnovaComponent, double> sgpp::base::OperationAnova::calculateAnovaOrderVariances(
    const DataVector& alpha) {
  return AnovaHelper::calculateAnovaOrderVariances(gridStorage, alpha);
}
