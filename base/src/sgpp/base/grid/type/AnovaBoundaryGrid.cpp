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

void AnovaBoundaryGrid::serialize(std::ostream& ostr, int version) {
  this->Grid::serialize(ostr, version);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& AnovaBoundaryGrid::getGenerator() { return generator; }

}  // namespace base
}  // namespace sgpp
