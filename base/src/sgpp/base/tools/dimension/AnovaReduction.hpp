#include <sgpp/base/tools/Sample.hpp>

namespace sgpp {
namespace base { 
namespace AnovaReduction {

  struct Result
  {
  SGridSample reducedSample;
    double coveredVariance;
  std::list<size_t> dimensionOrder;
    size_t activeDimensions;
  std::vector<double> sensitivities;
  };
  
    bool reduce(const SGridSample& sample, double coveredVariance, size_t maxEffectiveDims, bool removeInactiveDimensions, Result& result);
}
}
}  // namespace sgpp