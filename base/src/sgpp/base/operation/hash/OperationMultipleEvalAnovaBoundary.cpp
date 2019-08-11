// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalAnovaBoundary.hpp>
#include "common/basis/AnovaLinearBoundaryBasis.hpp"

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

void OperationMultipleEvalAnovaBoundary::mult(DataVector& alpha, DataVector& result) {
  AlgorithmDGEMV<SAnovaLinearBoundaryBasis> op;
  AnovaLinearBoundaryBasis<unsigned int, unsigned int> base;

  op.mult(storage, base, alpha, this->dataset, result);
}

void OperationMultipleEvalAnovaBoundary::multTranspose(DataVector& source, DataVector& result) {
  AlgorithmDGEMV<SAnovaLinearBoundaryBasis> op;
  AnovaLinearBoundaryBasis<unsigned int, unsigned int> base;

  op.mult_transposed(storage, base, source, this->dataset, result);
}

double OperationMultipleEvalAnovaBoundary::getDuration() { return 0.0; }

}  // namespace base
}  // namespace sgpp
