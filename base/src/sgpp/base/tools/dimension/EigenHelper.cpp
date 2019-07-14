#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include "EigenHelper.hpp"

sgpp::base::DataVector EigenHelper::fromEigen(const Eigen::VectorXd& e) {
  Eigen::VectorXd copy = e;
  sgpp::base::DataVector v(copy.data(), e.size());
  return std::move(v);
}

sgpp::base::DataMatrix EigenHelper::fromEigen(const Eigen::MatrixXd& e) {
  sgpp::base::DataMatrix m(e.data(), e.cols(), e.rows());
  m.transpose();
  return std::move(m);
}

Eigen::MatrixXd EigenHelper::toEigen(const sgpp::base::DataMatrix& matrix) {
  sgpp::base::DataMatrix copy = matrix;
  Eigen::Map<Eigen::MatrixXd> m(copy.getPointer(), matrix.getNcols(), matrix.getNrows());
  Eigen::MatrixXd mat = m.matrix();
  mat.transposeInPlace();
  return std::move(mat);
}

Eigen::VectorXd EigenHelper::toEigen(const sgpp::base::DataVector& vector) {
  sgpp::base::DataVector copy = vector;
  Eigen::Map<Eigen::VectorXd> v(copy.data(), vector.size());
  return std::move(v);
}

sgpp::base::DataMatrix EigenHelper::mult(sgpp::base::DataMatrix& m1, sgpp::base::DataMatrix& m2) {
  Eigen::MatrixXd e1 = toEigen(m1);
  Eigen::MatrixXd e2 = toEigen(m2);
  Eigen::MatrixXd e3 = e1 * e2;
  return fromEigen(e3);
}

sgpp::base::DataVector EigenHelper::mult(sgpp::base::DataMatrix& m,
                                         const sgpp::base::DataVector& v) {
  Eigen::MatrixXd em = toEigen(m);
  Eigen::VectorXd ev = toEigen(v);
  Eigen::VectorXd r = em * ev;
  return fromEigen(r);
}

void EigenHelper::svd(const Eigen::MatrixXd& input, sgpp::base::DataMatrix& eigenVectorMatrix,
                      sgpp::base::DataVector& eigenValues) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(input);
  if (eigensolver.info() != Eigen::Success) abort();
  const Eigen::MatrixXd& m = eigensolver.eigenvectors();
  eigenVectorMatrix = fromEigen(m);
  const Eigen::VectorXd& v = eigensolver.eigenvalues();
  eigenValues = fromEigen(v);

  size_t dimensions = input.cols();
  for (size_t c = 0; c < dimensions; c++) {
    size_t max = c;
    double maxValue = 0.0;
    for (size_t j = c; j < dimensions; j++) {
      if (eigenValues[j] > maxValue) {
        max = j;
        maxValue = eigenValues[j];
      }
    }
    sgpp::base::DataVector temp1(dimensions);
    eigenVectorMatrix.getColumn(c, temp1);
    sgpp::base::DataVector temp2(dimensions);
    eigenVectorMatrix.getColumn(max, temp2);
    eigenVectorMatrix.setColumn(c, temp2);
    eigenVectorMatrix.setColumn(max, temp1);
    double temp3 = eigenValues[c];
    eigenValues[c] = eigenValues[max];
    eigenValues[max] = temp3;
  }
}
