#ifdef USE_EIGEN

#include <eigen3/Eigen/Eigenvalues>
#include <sgpp/base/tools/EigenHelper.hpp>

namespace sgpp {
namespace base {

DataVector EigenHelper::fromEigen(const Eigen::VectorXd& e) {
  Eigen::VectorXd copy = e;
  DataVector v(copy.data(), e.size());
  return std::move(v);
}

DataMatrix EigenHelper::fromEigen(const Eigen::MatrixXd& e) {
  // Eigen uses column-major order, so we have to transpose
  DataMatrix m(e.data(), e.cols(), e.rows());
  m.transpose();
  return std::move(m);
}

Eigen::MatrixXd EigenHelper::toEigen(const DataMatrix& matrix) {
  // Eigen uses column-major order, so we have to transpose
  DataMatrix copy = matrix;
  Eigen::Map<Eigen::MatrixXd> m(copy.getPointer(), matrix.getNcols(), matrix.getNrows());
  Eigen::MatrixXd mat = m.matrix();
  mat.transposeInPlace();
  return std::move(mat);
}

Eigen::VectorXd EigenHelper::toEigen(const DataVector& vector) {
  DataVector copy = vector;
  Eigen::Map<Eigen::VectorXd> v(copy.data(), vector.size());
  return std::move(v);
}

DataMatrix EigenHelper::mult(const DataMatrix& m1, const DataMatrix& m2) {
  Eigen::MatrixXd e1 = toEigen(m1);
  Eigen::MatrixXd e2 = toEigen(m2);
  Eigen::MatrixXd e3 = e1 * e2;
  return fromEigen(e3);
}

DataVector EigenHelper::mult(const DataMatrix& m, const DataVector& v) {
  Eigen::MatrixXd em = toEigen(m);
  Eigen::VectorXd ev = toEigen(v);
  Eigen::VectorXd r = em * ev;
  return fromEigen(r);
}

void EigenHelper::svd(const Eigen::MatrixXd& input, DataMatrix& eigenVectorMatrix,
                      DataVector& eigenValues) {
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
    DataVector temp1(dimensions);
    eigenVectorMatrix.getColumn(c, temp1);
    DataVector temp2(dimensions);
    eigenVectorMatrix.getColumn(max, temp2);
    eigenVectorMatrix.setColumn(c, temp2);
    eigenVectorMatrix.setColumn(max, temp1);
    double temp3 = eigenValues[c];
    eigenValues[c] = eigenValues[max];
    eigenValues[max] = temp3;
  }
}


void EigenHelper::solveSLE(const sgpp::base::DataMatrix& A, const sgpp::base::DataVector& b,
  DataVector& x) {
  auto eigenA = toEigen(A);
  auto eigenb = toEigen(b);
  Eigen::VectorXd eigenx = eigenA.householderQr().solve(eigenb);
  x = fromEigen(eigenx);
}
}  // namespace base
}  // namespace sgpp

#endif
