// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationWeightedSecondMoment.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/globaldef.hpp>

#include <functional>
#include <map>
#include <tuple>

namespace sgpp {
namespace base {

/**
 * This class provides the second moment of a sparse grid function f given as interpolant of any not
 * a knot B-spline basis w.r.t. a probability density function \rho, i.e. \int
 * f(x)^2 \rho (x) dx
 */
class OperationWeightedSecondMomentNak : public OperationWeightedSecondMoment {
 public:
  /**
   * Constructor of OperationQuadratureNakBsplineExtended
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   */
  OperationWeightedSecondMomentNak(GridStorage& storage, GridType gridType, size_t degree)
      : storage(storage), degree(degree) {
    basis = initializeBasis(gridType, degree);
  }

  ~OperationWeightedSecondMomentNak() override {}

  /**
   * initializes a nak B spline basis according to gridType and degree
   * @param gridType	type of the basis
   * @param degree		degree of the basis
   * @return 			the nak Bsplien basis
   */
  std::shared_ptr<sgpp::base::SBasis> initializeBasis(sgpp::base::GridType gridType, size_t degree);

  /**
   * Quadrature for not a knot B-spline basis functions w.r.t. a probability density function
   *
   * @param alpha   	Coefficient vector for current grid
   * @param pdf			probability density function
   * @parm quadOrder	order for the gauss Legendre quadrature
   */
  double doWeightedQuadrature(DataVector& alpha, std::shared_ptr<sgpp::base::Distribution> pdf,
                              size_t quadOrder);

 private:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  // degree of the nak B-spline basis
  size_t degree;
  /// NakBsplineBoundaryCombigrid Basis object
  std::shared_ptr<sgpp::base::SBasis> basis;
  // quadrature coordinates
  sgpp::base::DataVector coordinates;
  // quadrature weights
  sgpp::base::DataVector weights;
  // tuple used as hash to store scalar products in innerProducts
  typedef std::tuple<size_t, size_t, size_t, size_t> hashType;
  // hash storage for scalar products. Holds all calculated scalar products s.t. they do not have to
  // calculated again if the same combination of indices and levels  is queried
  std::map<hashType, double> innerProducts;

  /**
   * calculates the one dimensional integral \int f*g \rho dx where f and g are B-spline basis
   * functions  and \rho is a probability  DensityFunction
   *
   * @param level1 	level of the first B-spline
   * @param index1 	index of the first B-spline
   * @param level2 	level of the second B-spline
   * @param index2 	index of the second B-spline
   * @param pdf		probability density function
   *
   * @return  integral first basis function *  second basis function * probability density function
   */
  double weightedBasisScalarProduct(unsigned int level1, unsigned int index1, unsigned int level2,
                                    unsigned int index2,
                                    std::shared_ptr<sgpp::base::Distribution> pdf);

  /**
   * used to get the support segments of a not a knot B-spline basis functions.
   *
   * @param level	level of the B-spline basis function
   * @param index	index of the B-spline basis function
   * @param degree	degree of the B-spline basis function
   *
   * @return the indices of the segments of the not a knot B-spline basis functions support
   */
  sgpp::base::DataVector nakBSplineSupport(size_t level, size_t index, size_t degree);
};
}  // namespace base
}  // namespace sgpp
