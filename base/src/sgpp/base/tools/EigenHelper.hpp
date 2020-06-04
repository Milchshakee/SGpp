// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#ifdef USE_EIGEN

#include <eigen3/Eigen/Core>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp {
namespace base {

/**
 * Provides useful functionality for vectors and matrices by using functions of the Eigen library.
 */
namespace EigenHelper {
/**
 * Converts an Eigen vector object to an SGpp vector object.
 */
sgpp::base::DataVector fromEigen(const Eigen::VectorXd& e);

/**
 * Converts an Eigen matrix object to an SGpp matrix object.
 */
sgpp::base::DataMatrix fromEigen(const Eigen::MatrixXd& e);

/**
 * Converts an SGpp vector object to an Eigen vector object.
 */
Eigen::MatrixXd toEigen(const sgpp::base::DataMatrix& matrix);

/**
 * Converts an SGpp matrix object to an Eigen matrix object.
 */
Eigen::VectorXd toEigen(const sgpp::base::DataVector& vector);

/**
 * Multiplies two matrices.
 */
sgpp::base::DataMatrix mult(const sgpp::base::DataMatrix& m1, const sgpp::base::DataMatrix& m2);

/**
 * Multiplies a matrix with a vector.
 */
sgpp::base::DataVector mult(const sgpp::base::DataMatrix& m, const sgpp::base::DataVector& v);

/**
 * Performs a singular value decomposition on a matrix.
 */
void svd(const Eigen::MatrixXd& input, sgpp::base::DataMatrix& eigenVectors,
         sgpp::base::DataVector& eigenValues);

  void solveSLE(const sgpp::base::DataMatrix& A, const sgpp::base::DataVector& b, DataVector& x);
}  // namespace EigenHelper
}  // namespace base
}  // namespace sgpp

#endif
