// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/type/AnovaBoundaryGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/AnovaBoundaryBasis.hpp>
#include <ostream>

namespace sgpp {
namespace base {


AnovaBoundaryGrid::AnovaBoundaryGrid(std::istream& istr)
    : Grid(istr), generator(storage) 
  {}

AnovaBoundaryGrid::AnovaBoundaryGrid(size_t dim)
    : Grid(dim), generator(storage) {}

sgpp::base::GridType AnovaBoundaryGrid::getType() { return sgpp::base::GridType::AnovaBoundary; }

SBasis& AnovaBoundaryGrid::getBasis() {
  static SAnovaBoundaryBasis basis;
  return basis;
}

Grid* AnovaBoundaryGrid::unserialize(std::istream& istr)
{
  return new AnovaBoundaryGrid(istr);
}

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
