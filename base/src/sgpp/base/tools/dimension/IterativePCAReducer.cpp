#include <sgpp/base/tools/dimension/IterativePCAReducer.hpp>
#include "Tools.hpp"


  std::unique_ptr<sgpp::optimization::ScalarFunction> IterativePCAReducer::reduceData(
    sgpp::base::DataMatrix& input) {
  size_t dimension = input.getNcols();
      //Tools::centerMean(x);
  sgpp::base::DataMatrix random(dimension, dimension);
  for (size_t d = 0; d < dimension; ++d) {
    sgpp::base::DataVector rand;
    random.setColumn(d, rand);
    }

  Eigen::MatrixXd x = Tools::toEigen(input);
  Eigen::MatrixXd r = Tools::toEigen(random);
  Eigen::MatrixXd s(dimension, dimension);
  Eigen::MatrixXd e(dimension, dimension);
  for (size_t it = 0; it < iterations; ++it) {
    s = x.transpose() * (x * r);
    e = r.transpose() * s;
    r = s;
  }

    sgpp::base::DataMatrix eigenVectorMatrix(dimension, dimension);
  sgpp::base::DataVector eigenValues(dimension);
  for (size_t c = 0; c < dimension; c++) {
    eigenValues.set(c, e(c, c));
  }
  return nullptr;
  }

