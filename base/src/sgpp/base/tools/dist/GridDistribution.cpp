#include <sgpp/base/tools/dist/GridDistribution.hpp>

namespace sgpp {
namespace base {

GridDistribution::GridDistribution(sgpp::base::Grid& grid)
    : VectorDistribution(grid.getSize(), grid.getDimension()) {
  for (size_t i = 0; i < grid.getSize(); i++) {
    sgpp::base::GridPoint& gp = grid.getStorage().getPoint(i);
    for (size_t j = 0; j < grid.getDimension(); j++) {
      vectors[i][j] = gp.getStandardCoordinate(j);
    }
  }
}
}  // namespace base
}  // namespace sgpp
