#include <sgpp/base/operation/hash/OperationVariance.hpp>
#include "sgpp/base/grid/GridStorage.hpp"
#include "common/basis/AnovaBoundaryBasis.hpp"

double getIntegral(sgpp::base::Grid& grid, const sgpp::base::OperationVariance::LevelVector& levels) {
  sgpp::base::SAnovaBoundaryBasis& b = dynamic_cast<sgpp::base::SAnovaBoundaryBasis&>(grid.getBasis());
  double val = 1;
  for (size_t d = 0; d < levels.size(); d++) {
    val *= b.getIntegral(levels[d]);
  }
  return val;
}

void sgpp::base::OperationVariance::calculateIncrementVariance(sgpp::base::Grid& grid,
                                                               const DataVector& alpha,
                                                               const LevelVector& levels) {
  const sgpp::base::GridStorage& gridStorage = grid.getStorage();
  out: for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      if (gp.getLevel(d) != levels[d]) {
        goto out;
      }
    }
    double val = alpha[i];
    double 
  }

  for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
    result += iter->second * alpha[iter->first];
  }

  return result;
}
