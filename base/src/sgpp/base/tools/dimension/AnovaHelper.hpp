#ifndef ANOVAHELPER
#define ANOVAHELPER

#include <vector>
#include "sgpp/base/grid/storage/hashmap/HashGridPoint.hpp"
#include "sgpp/base/tools/Sample.hpp"

namespace sgpp {
namespace base {
namespace AnovaHelper {

/**
 * Vector that holds levels for every dimension
 */
typedef std::vector<sgpp::base::HashGridPoint::level_type> LevelVector;

typedef std::vector<bool> AnovaComponent;
typedef std::vector<AnovaComponent> AnovaComponentVector;

AnovaComponent getAnovaComponentOfPoint(const GridPoint& point);

Sample<AnovaComponent, double> calculateAnovaOrderVariances(GridStorage& gridStorage,
                                                            const DataVector& alpha);

}  // namespace AnovaHelper
}  // namespace base
}  // namespace sgpp

#endif