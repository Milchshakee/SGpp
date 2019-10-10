// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <ostream>
#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

namespace sgpp {
namespace base {

AnovaBoundaryGrid::AnovaBoundaryGrid(std::istream& istr) : Grid(istr), generator(storage) {}

AnovaBoundaryGrid::AnovaBoundaryGrid(size_t dim) : Grid(dim), generator(storage) {}

AnovaBoundaryGrid::AnovaBoundaryGrid(size_t dim, std::vector<AnovaTypes::LevelIndexPair>& anchor)
    : Grid(dim), generator(storage), anchor(anchor) {}

void AnovaBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& AnovaBoundaryGrid::getGenerator() { return generator; }

bool AnovaBoundaryGrid::hasAnchor() { return anchor.size() > 0; }

const std::vector<AnovaTypes::LevelIndexPair>& AnovaBoundaryGrid::getAnchor() { return anchor; }

std::shared_ptr<ScalarFunction> AnovaBoundaryGrid::getSamplingFunction(ScalarFunction& func) {
  if (!hasAnchor()) {
    throw operation_exception("Grid has no anchor"); 
    }

  std::function<double(const DataVector&)> f = [&func, this](const DataVector& v) {
    DataVector out = v;
    for (size_t d = 0; d < getDimension(); d++) {
      if (v[d] == 0 && anchor[d].level > -1) {
        out[d] = static_cast<double>(anchor[d].index) / static_cast<double>(1 << anchor[d].level);
      }
    }
    return func.eval(out);
  };

  return std::make_shared<WrapperScalarFunction>(func.getNumberOfParameters(), f);
}
}  // namespace base
}  // namespace sgpp
