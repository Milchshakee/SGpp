// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>
#include <vector>
#include <sgpp/base/operation/hash/OperationEvalAnovaLinearBoundary.hpp>

namespace sgpp {
namespace base {

double OperationEvalAnovaLinearBoundary::eval(const DataVector& alpha, const DataVector& point) {
  typedef std::vector<std::pair<size_t, double> > IndexValVector;

  IndexValVector vec;
  AnovaLinearBoundaryBasis<unsigned int, unsigned int> base;
  GetAffectedBasisFunctions<AnovaLinearBoundaryBasis<unsigned int, unsigned int> > ga(storage);

  ga(base, point, vec);

  double result = 0.0;

  for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
    GridPoint& p = storage.getPoint(iter->first);
    if (component.empty() || AnovaBoundaryGrid::getAnovaComponentOfPoint(p) == component) {
      result += iter->second * alpha[iter->first];
    }
  }

  return result;
}

}  // namespace base
}  // namespace sgpp
