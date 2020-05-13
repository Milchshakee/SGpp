#pragma once

#include <sgpp/base/tools/Sample.hpp>


namespace sgpp {
namespace base {
namespace AnovaAnalyser
{

sgpp::base::SGridSample createReducedAnovaSample(sgpp::base::SGridSample& sample, size_t reducedDims);

struct Info {
  std::list<size_t> dimensionOrder;
  std::vector<double> totalEffectIndices;
};

Info evaluateFunction(sgpp::base::SGridSample& sample, DistributionSample& dist);

  size_t getReducedDimensions(Info& i, double minVariance);

  class TransformationFunction : public VectorFunction
  {

  public:

    TransformationFunction(size_t d, const std::vector<size_t>& dimension_order,
                          size_t reduced_dims);
    ~TransformationFunction() override;
    void eval(const DataVector& x, DataVector& value) override;
    void clone(std::unique_ptr<VectorFunction>& clone) const override;

  private:
    std::vector<size_t> dimensionOrder;
  };
  }
}  // namespace base
}  // namespace sgpp