// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

BoundingBox::BoundingBox(size_t dimension) :
    dimension(dimension),
    boundingBox1Ds(dimension, {0.0, 1.0}) {
}

BoundingBox::BoundingBox(const std::vector<BoundingBox1D>& boundingBox1Ds) :
    dimension(boundingBox1Ds.size()),
    boundingBox1Ds(boundingBox1Ds) {
}

void BoundingBox::setBoundary(size_t d, const BoundingBox1D& boundingBox1D) {
  boundingBox1Ds[d] = boundingBox1D;
}

BoundingBox1D BoundingBox::getBoundary(size_t d) const {
  return boundingBox1Ds[d];
}

size_t BoundingBox::getDimensions() const {
  return dimension;
}

double BoundingBox::getIntervalWidth(size_t d) const {
  return boundingBox1Ds[d].rightBoundary - boundingBox1Ds[d].leftBoundary;
}

double BoundingBox::getIntervalOffset(size_t d) const {
  return boundingBox1Ds[d].leftBoundary;
}

bool BoundingBox::isUnitCube() const {
  for (size_t d = 0; d < dimension; d++) {
    if ((boundingBox1Ds[d].leftBoundary != 0.0) ||
        (boundingBox1Ds[d].rightBoundary != 1.0)) {
      return false;
    }
  }

  return true;
}

bool BoundingBox::hasDirichletBoundaryLeft(size_t d) const {
  return boundingBox1Ds[d].bDirichletLeft;
}

bool BoundingBox::hasDirichletBoundaryRight(size_t d) const {
  return boundingBox1Ds[d].bDirichletRight;
}

void BoundingBox::toString(std::string& text) const {
  std::stringstream str;

  for (size_t d = 0; d < dimension; d++) {
    str << "Dimensions " << d << "(" << boundingBox1Ds[d].leftBoundary
        << "," <<  boundingBox1Ds[d].rightBoundary << ")\n";
  }

  text = str.str();
}

std::string BoundingBox::toString() const {
  std::string str;
  toString(str);
  return str;
}

}  // namespace base
}  // namespace sgpp
