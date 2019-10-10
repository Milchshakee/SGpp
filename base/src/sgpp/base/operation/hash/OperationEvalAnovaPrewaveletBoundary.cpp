#include <sgpp/base/operation/hash/OperationEvalAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/operation/hash/common/basis/AnovaPrewaveletBoundaryBasis.hpp>

double sgpp::base::OperationEvalAnovaPrewaveletBoundary::eval(const DataVector& alpha,
                                                              const DataVector& point)
{
  typedef std::vector<std::pair<size_t, double>> IndexValVector;

  IndexValVector vec;
  AnovaGridIterator it(storage);

  static SAnovaPrewaveletBoundaryBasis basis;

  for (size_t i = 0; i < storage.getSize(); i++) {
    GridPoint& p = storage.getPoint(i);
    if (component.empty() || AnovaBoundaryGrid::getAnovaComponentOfPoint(p) == component) {
      double value = 1.0;
      for (size_t d = 0; d < storage.getDimension(); ++d) {
        index_t current_index;
        sgpp::base::AnovaBoundaryGrid::level_t current_level;
        AnovaBoundaryGrid::fromNormalGridPointLevelIndex(p.getLevel(d), p.getIndex(d),
                                                         current_level, current_index);
        value *= current_level == -1 ? 1 : basis.eval(current_level, current_index, point[d]);
      }
      vec.push_back(std::make_pair(storage.getSequenceNumber(p), value));
    }
  }

  double result = 0.0;
  for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
    result += iter->second * alpha[iter->first];
  }

  return result;
}
