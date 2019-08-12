#include <sgpp/base/operation/hash/OperationEvalAnovaPrewaveletBoundary.hpp>
#include <sgpp/base/operation/hash/common/basis/AnovaPrewaveletBoundaryBasis.hpp>

void sgpp::base::OperationEvalAnovaPrewaveletBoundary::rec(
    const sgpp::base::DataVector& point, size_t current_dim, sgpp::base::AnovaGridIterator& iter,
    std::vector<std::pair<size_t, double>>& result) {
  static SAnovaPrewaveletBoundaryBasis basis;
  double value = 1.0;
  for (size_t d = 0; d < storage.getDimension(); ++d) {
    index_t current_index;
    sgpp::base::AnovaBoundaryGrid::level_t current_level;
    iter.get(d, current_level, current_index);
    value *= current_level == -1 ? 1 : basis.eval(current_level, current_index, point[d]);
  }
  result.push_back(std::make_pair(iter.seq(), value));

  for (size_t d = current_dim; d < storage.getDimension(); d++) {
    index_t current_index;
    sgpp::base::AnovaBoundaryGrid::level_t current_level;
    iter.get(d, current_level, current_index);

    if (!component.empty() && !component[d]) {
      //Do nothing
    } else if (current_level == -1) {
      iter.resetToLevelZeroInDim(d);
      if (!storage.isInvalidSequenceNumber(iter.seq())) {
        rec(point, d, iter, result);
      }
    } else if (current_level == 0) {
      iter.resetToLevelOneInDim(d);
      if (!storage.isInvalidSequenceNumber(iter.seq())) {
        rec(point, d, iter, result);
      }
    } else {
      if (iter.hintLeft(d)) {
        iter.leftChild(d);
        rec(point, d, iter, result);
        iter.up(d);
      }

      if (iter.hintRight(d)) {
        iter.rightChild(d);
        rec(point, d, iter, result);
        iter.up(d);
      }
    }
  }
}

double sgpp::base::OperationEvalAnovaPrewaveletBoundary::eval(const DataVector& alpha,
                                                              const DataVector& point)
{
  typedef std::vector<std::pair<size_t, double>> IndexValVector;

  IndexValVector vec;
  AnovaGridIterator it(storage);
  //No need to calculate affected basis functions since every function is affected
  rec(point, 0, it, vec);

  double result = 0.0;
  for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
    result += iter->second * alpha[iter->first];
  }

  return result;
}
