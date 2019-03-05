// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/activeSubspaces/SparseGridResponseSurfaceNakBspline.hpp>

namespace sgpp {
namespace datadriven {

void SparseGridResponseSurfaceNakBspline::initialize() {
  numDim = objectiveFunc->getNumberOfParameters();
  if (gridType == sgpp::base::GridType::NakBspline) {
    grid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    grid = std::make_shared<sgpp::base::NakBsplineModifiedGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
    grid = std::make_shared<sgpp::base::NakBsplineExtendedGrid>(numDim, degree);
    basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
  } else {
    throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
  }
}

void SparseGridResponseSurfaceNakBspline::regular(size_t level) {
  grid->getGenerator().regular(level);
  calculateInterpolationCoefficients();
  interpolant =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

void SparseGridResponseSurfaceNakBspline::regularByPoints(size_t numPoints) {
  size_t level = 1;
  do {
    // todo (rehmemk) instead of trying out until pointnumber matches, use formula for number of
    // grid points
    grid->getStorage().clear();
    grid->getGenerator().regular(level);
    level++;
  } while (grid->getSize() < numPoints);
  calculateInterpolationCoefficients();
  interpolant =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

void SparseGridResponseSurfaceNakBspline::surplusAdaptive(size_t maxNumGridPoints,
                                                          size_t initialLevel,
                                                          size_t refinementsNum) {
  regular(initialLevel);
  while (grid->getSize() < maxNumGridPoints) {
    refineSurplusAdaptive(refinementsNum);
    interpolant =
        std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
    interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
        *grid, coefficients);
  }
}

void SparseGridResponseSurfaceNakBspline::regularData(size_t level,
                                                      sgpp::base::DataMatrix evaluationPoints,
                                                      sgpp::base::DataVector functionValues,
                                                      double lambda) {
  grid->getGenerator().regular(level);
  double mse = 0;
  sgpp::base::DataVector errorPerBasis;
  coefficients = EigenRegression(grid, degree, DataMatrixToEigen(evaluationPoints), functionValues,
                                 mse, errorPerBasis);
  interpolant =
      std::make_unique<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_unique<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

double SparseGridResponseSurfaceNakBspline::eval(sgpp::base::DataVector v) {
  return interpolant->eval(v);
}
double SparseGridResponseSurfaceNakBspline::evalGradient(sgpp::base::DataVector v,
                                                         sgpp::base::DataVector& gradient) {
  return interpolantGradient->eval(v, gradient);
}

double SparseGridResponseSurfaceNakBspline::evalNonUniform(sgpp::base::DataVector v,
                                                           sgpp::base::DataVector lBounds,
                                                           sgpp::base::DataVector uBounds) {
  sgpp::base::DataVector newlBounds(lBounds.getSize(), 0.0);
  sgpp::base::DataVector newuBounds(lBounds.getSize(), 1.0);
  transformPoint(v, lBounds, uBounds, newlBounds, newuBounds);
  return interpolant->eval(v);
}

double SparseGridResponseSurfaceNakBspline::evalGradientNonUniform(sgpp::base::DataVector v,
                                                                   sgpp::base::DataVector& gradient,
                                                                   sgpp::base::DataVector lBounds,
                                                                   sgpp::base::DataVector uBounds) {
  sgpp::base::DataVector newlBounds(lBounds.getSize(), 0.0);
  sgpp::base::DataVector newuBounds(lBounds.getSize(), 1.0);
  transformPoint(v, lBounds, uBounds, newlBounds, newuBounds);
  return interpolantGradient->eval(v, gradient);
}

double SparseGridResponseSurfaceNakBspline::getIntegral() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  double integral = 0;
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    double integral1D = 1;
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      integral1D *=
          basis->getIntegral(gridStorage.getPointLevel(i, d), gridStorage.getPointIndex(i, d));
    }
    integral += coefficients[i] * integral1D;
  }
  return integral;
}

// ----------------- auxiliary routines -----------

void SparseGridResponseSurfaceNakBspline::refineSurplusAdaptive(size_t refinementsNum) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  calculateInterpolationCoefficients();
}

void SparseGridResponseSurfaceNakBspline::calculateInterpolationCoefficients() {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector f_values(gridStorage.getSize(), 0.0);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      p[d] = gridStorage.getPointCoordinate(i, d);
    }
    f_values[i] = objectiveFunc->eval(p);
  }
  sgpp::optimization::sle_solver::Auto sleSolver;
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  if (!sleSolver.solve(hierSLE, f_values, coefficients)) {
    std::cout << "Solving failed!" << std::endl;
  }
}

void SparseGridResponseSurfaceNakBspline::transformPoint(sgpp::base::DataVector& v,
                                                         sgpp::base::DataVector lBounds,
                                                         sgpp::base::DataVector uBounds,
                                                         sgpp::base::DataVector newlBounds,
                                                         sgpp::base::DataVector newuBounds) {
  v.sub(lBounds);
  uBounds.sub(lBounds);
  v.componentwise_div(uBounds);
  newuBounds.sub(newlBounds);
  v.componentwise_mult(newuBounds);
  v.add(newlBounds);
}

}  // namespace datadriven
}  // namespace sgpp
