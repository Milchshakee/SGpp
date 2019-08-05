#include <sgpp/base/tools/dist/FixedDistribution.hpp>

namespace sgpp {
namespace base {

FixedDistribution::FixedDistribution(const DataMatrix& data)
    : VectorDistribution(data.getNrows(), data.getNcols()) {
  for (size_t i = 0; i < size; i++) {
    data.getRow(i, vectors[i]);
  }
}

FixedDistribution::FixedDistribution(size_t size, size_t dimensions,
                                     const std::vector<DataVector>& vectors)
    : VectorDistribution(size, dimensions) {
  VectorDistribution::vectors = vectors;
}
}  // namespace base
}  // namespace sgpp