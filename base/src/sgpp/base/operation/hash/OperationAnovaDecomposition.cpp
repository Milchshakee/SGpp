#include <sgpp/base/algorithm/sweep_anova.hpp>
#include <sgpp/base/grid/storage/hashmap/AnovaGridIterator.hpp>
#include <sgpp/base/operation/hash/OperationAnovaDecomposition.hpp>

size_t getIndex(sgpp::base::AnovaGridIterator& it,
                const sgpp::base::AnovaBoundaryGrid::AnovaComponent& component) {
  sgpp::base::AnovaGridIterator copy = it;
  for (size_t dim = 0; dim < component.size(); dim++) {
    if (!component[dim]) {
      copy.set(dim, -1, 0);
    }
  }
  return copy.seq();
}

void sgpp::base::OperationAnovaDecomposition::decompose(const DataVector& alpha,
                                                        DataVector& result) {}

void sgpp::base::OperationAnovaDecomposition::calcExpectedValue(
    const DataVector& alpha, DataVector& result,
    const AnovaBoundaryGrid::AnovaComponent& component) {
  typedef std::function<void(DataVector&, DataVector&, AnovaGridIterator&,
                             const AnovaBoundaryGrid::AnovaComponent&)>
      FuncType;

  FuncType f = [this, &component, &alpha](DataVector& in, DataVector& out,
                                                   AnovaGridIterator& it,
                                                   const AnovaBoundaryGrid::AnovaComponent& c) {
    double funcMean = 1;
    for (size_t dim = 0; dim < storage.getDimension(); dim++) {
      if (component[dim]) {
        AnovaTypes::level_t level;
        AnovaTypes::index_t index;
        it.get(dim, level, index);
        double mean;

        if (level == -1) {
          mean = 1. / 2.;
        } else if (level == 0) {
          mean = 1. / 6.;
        } else {
          // evaluate the first moment in the unit hypercube
          level_t l;
          index_t i;
          AnovaBoundaryGrid::toNormalGridPointLevelIndex(level, index, l, i);
          mean = static_cast<double>(i) * std::pow(4.0, -static_cast<double>(l));
        }
        funcMean *= mean;
      }
    }

    out[getIndex(it, component)] += alpha.get(it.seq()) * funcMean;
  };
  sweep_anova_component<FuncType> sweep(f, storage);
  sweep.sweep1D_AnovaBoundary(alpha, result, component);
}
