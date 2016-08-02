// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPOINTORDERING_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPOINTORDERING_HPP_

#include <sgpp/combigrid/grid/ordering/AbstractPointOrdering.hpp>

namespace sgpp {
namespace combigrid {

class ExponentialLevelorderPointOrdering : public AbstractPointOrdering {
 public:
  virtual ~ExponentialLevelorderPointOrdering();

  virtual size_t convertIndex(size_t level, size_t numPoints, size_t index);

  virtual size_t numPoints(size_t level);

  virtual std::shared_ptr<AbstractPermutationIterator> getSortedPermutationIterator(
      size_t level, std::vector<double> const &points, size_t numPoints);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPOINTORDERING_HPP_ \
          */