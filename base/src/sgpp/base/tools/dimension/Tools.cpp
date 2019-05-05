#include "Tools.hpp"
#include <gsl/gsl_eigen.h>

sgpp::base::DataMatrix Tools::fromEigen(Eigen::MatrixXd& e) {
  sgpp::base::DataMatrix m(e.data(), e.rows(), e.cols());
  return std::move(m);
}

Eigen::MatrixXd Tools::toEigen(sgpp::base::DataMatrix& matrix) {
  Eigen::Map<Eigen::MatrixXd> m(matrix.getPointer(), matrix.getNrows(), matrix.getNcols());
  return std::move(m.matrix());
}

sgpp::base::DataMatrix Tools::mult(sgpp::base::DataMatrix& m1,
                                   sgpp::base::DataMatrix& m2) {
  Eigen::MatrixXd e1 = toEigen(m1);
  Eigen::MatrixXd e2 = toEigen(m2);
  Eigen::MatrixXd e3 = e1 * e2;
  return fromEigen(e3);
}

sgpp::base::DataMatrix Tools::centerMean(sgpp::base::DataMatrix& matrix) {
}


void Tools::svd(const sgpp::base::DataMatrix& input, sgpp::base::DataMatrix& eigenVectorMatrix,
  sgpp::base::DataVector& eigenValues) {
  size_t dimensions = input.getNcols();
  sgpp::base::DataMatrix copy = input;
  gsl_matrix_view m = gsl_matrix_view_array(copy.getPointer(), dimensions, dimensions);
  auto q = std::unique_ptr<gsl_matrix>{gsl_matrix_alloc(dimensions, dimensions)};
  auto e = std::unique_ptr<gsl_vector>{gsl_vector_alloc(dimensions)};
  gsl_eigen_symmv_workspace* ws = gsl_eigen_symmv_alloc(dimensions);
  gsl_eigen_symmv(&m.matrix, e.get(), q.get(), ws);
  gsl_eigen_symmv_free(ws);

  for (size_t r = 0; r < dimensions; r++) {
    for (size_t c = 0; c < dimensions; c++) {
      eigenVectorMatrix.set(r, c, gsl_matrix_get(q.get(), r, c));
    }
  }
  for (size_t c = 0; c < dimensions; c++) {
    eigenValues.set(c, gsl_vector_get(e.get(), c));
  }

  for (size_t c = 0; c < dimensions; c++) {
    size_t max = 0;
    double maxValue = 0.0;
    for (size_t j = c; j < dimensions; j++) {
      if (eigenValues[j] > maxValue) {
        max = j;
        maxValue = eigenValues[j];
      }
    }
    sgpp::base::DataVector temp1;
    eigenVectorMatrix.getRow(c, temp1);
    sgpp::base::DataVector temp2;
    eigenVectorMatrix.getRow(max, temp2);
    eigenVectorMatrix.setColumn(c, temp2);
    eigenVectorMatrix.setColumn(max, temp1);
    double temp3 = eigenValues[c];
    eigenValues[c] = eigenValues[max];
    eigenValues[max] = temp3;
  }
  }
