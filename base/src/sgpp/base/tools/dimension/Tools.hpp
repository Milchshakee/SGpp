// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <eigen3/Eigen/Core>
#include "sgpp/base/datatypes/DataMatrix.hpp"

namespace Tools
{
sgpp::base::DataMatrix fromEigen(Eigen::MatrixXd& e);
Eigen::MatrixXd toEigen(sgpp::base::DataMatrix& matrix);
sgpp::base::DataMatrix mult(sgpp::base::DataMatrix& m1, sgpp::base::DataMatrix& m2);
sgpp::base::DataMatrix centerMean(sgpp::base::DataMatrix& matrix);
void svd(const sgpp::base::DataMatrix& input, sgpp::base::DataMatrix& eigenVectors,
         sgpp::base::DataVector& eigenValues);
}

#endif
