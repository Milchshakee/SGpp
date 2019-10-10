#include <sgpp/base/algorithm/sweep_anova.hpp>
#include <sgpp/base/grid/storage/hashmap/AnovaGridIterator.hpp>
#include <sgpp/base/operation/hash/OperationAnovaDecomposition.hpp>

using namespace sgpp::base;

typedef std::function<void(const sgpp::base::DataVector&, sgpp::base::DataVector&, sgpp::base::AnovaGridIterator&,
                           const sgpp::base::AnovaBoundaryGrid::AnovaComponent&)>
    FuncType;

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

void calcExpectedValue(const DataVector& alpha, DataVector& result, AnovaGridIterator& it,
                       const AnovaBoundaryGrid::AnovaComponent& component) {
  double sum = 0;

  FuncType f = [&sum, &component, &alpha](const DataVector& in, DataVector& out, AnovaGridIterator& it,
                                    const AnovaBoundaryGrid::AnovaComponent& c) {
    size_t levelSum = 0;
    size_t level0Count = 0;
    for (size_t dim = 0; dim < it.getStorage().getDimension(); dim++) {
      if (c[dim]) {
        AnovaTypes::level_t level;
        AnovaTypes::index_t index;
        it.get(dim, level, index);
        if (level == -1) {
          continue;
          }
        if (level == 0) {
          level0Count++;
        } else {
          levelSum += level;
          }
      }
    }

    double funcMean = pow(2.0, -static_cast<double>(levelSum));
    funcMean *= (pow(2.0, -static_cast<double>(level0Count)));

    //for (size_t dim = 0; dim < it.getStorage().getDimension(); dim++) {
    //  if (c[dim]) {
    //    AnovaTypes::level_t level;
    //    AnovaTypes::index_t index;
    //    it.get(dim, level, index);
    //    double mean;

    //    if (level == -1) {
    //      mean = 1. / 2.;
    //    } else if (level == 0) {
    //      mean = 1. / 6.;
    //    } else {
    //      // evaluate the first moment in the unit hypercube
    //      level_t l;
    //      index_t i;
    //      AnovaBoundaryGrid::toNormalGridPointLevelIndex(level, index, l, i);
    //      mean = static_cast<double>(i) * std::pow(4.0, -static_cast<double>(l));
    //    }
    //    funcMean *= mean;
    //  }
    //}
    std::cout << it.seq() << std::endl;
    std::cout << (funcMean) << std::endl;
    std::cout << (alpha.get(it.seq()) * funcMean) << std::endl;
    std::cout << getIndex(it, component) << std::endl;
    sum += alpha.get(it.seq()) * funcMean;
    //out[getIndex(it, component)] += alpha.get(it.seq()) * funcMean;
  };

  AnovaBoundaryGrid::AnovaComponent rest(component.size());
  for (size_t dim = 0; dim < component.size(); dim++) {
    rest[dim] = !component[dim];
  }

  sweep_anova_component<FuncType> sweep(f, it.getStorage());
  sweep.sweep1D_AnovaBoundary_Component(alpha, result, it, rest);

  
  double value = alpha[it.seq()];
  double change = sum - value;
  size_t currentSeq = it.seq();
  result[currentSeq] += change;
  
  // FuncType s = [currentSeq, change, &component, &alpha](const DataVector& in, DataVector& out, AnovaGridIterator& it,
  //                                  const AnovaBoundaryGrid::AnovaComponent& c) {
  //  if (currentSeq == it.seq()) {
  //    out[it.seq()] += change;
  //  } else {
  //    out[it.seq()] -= change;
  //  }
  //};

  //sweep_anova_component<FuncType> update(s, it.getStorage());
  //update.sweep1D_AnovaBoundary_Component(alpha, result, it, rest);
}

void sgpp::base::OperationAnovaDecomposition::decompose(DataVector& alpha,
                                                        DataVector& result) {
  FuncType func = calcExpectedValue;
  sweep_anova_component<FuncType> sweep(func, storage);
  sweep.sweep1D_AnovaBoundary_AllComponents(alpha, result);
}
