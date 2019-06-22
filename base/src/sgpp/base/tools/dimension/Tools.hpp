// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <eigen3/Eigen/Core>
#include "sgpp/base/datatypes/DataMatrix.hpp"

namespace Tools {
sgpp::base::DataVector fromEigen(const Eigen::VectorXd& e);
sgpp::base::DataMatrix fromEigen(const Eigen::MatrixXd& e);
Eigen::MatrixXd toEigen(const sgpp::base::DataMatrix& matrix);
Eigen::VectorXd toEigen(const sgpp::base::DataVector& vector);
sgpp::base::DataMatrix mult(sgpp::base::DataMatrix& m1, sgpp::base::DataMatrix& m2);
sgpp::base::DataVector mult(sgpp::base::DataMatrix& m, const sgpp::base::DataVector& v);
void svd(const sgpp::base::DataMatrix& input, sgpp::base::DataMatrix& eigenVectors,
         sgpp::base::DataVector& eigenValues);
}

#endif
