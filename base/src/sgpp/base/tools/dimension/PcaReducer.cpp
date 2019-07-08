#include <sgpp/base/tools/dimension/PcaReducer.hpp>
#include "Tools.hpp"


sgpp::base::PcaResult::PcaResult(const DataMatrix& m, size_t n, double coveredVariance) : coveredVariance(coveredVariance) {
  transformation = m;
  transformation.resizeRowsCols(n, m.getNcols());
}


sgpp::base::PcaFixedCutter::PcaFixedCutter(size_t n) : n(n) {
}

sgpp::base::PcaResult sgpp::base::PcaFixedCutter::cut(const VectorDistribution& input,
                                                      const PcaInfo& info){
  double sum = 0;
  for (size_t d = 0; d < n; ++d) {
    sum += info.varianceShares[d];
  }
  return PcaResult(info.basis, n, sum);
}

sgpp::base::PcaVarianceCutter::PcaVarianceCutter(double varianceShare) : minVarianceShare(varianceShare) {
}

sgpp::base::PcaResult sgpp::base::PcaVarianceCutter::cut(const VectorDistribution& input,
                                                         const PcaInfo& info) {
  double sum = 0;
  for (size_t d = 0; d < info.activeComponentsCount; ++d) {
    sum += info.varianceShares[d];
    if (sum >= minVarianceShare) {
      return PcaResult(info.basis, d + 1, sum);
    }
  }
  return PcaResult(info.basis, info.activeComponentsCount, sum);
}

sgpp::base::FixedDistribution sgpp::base::PcaResult::apply(const VectorDistribution& input) {
  std::vector<DataVector> newDist(input.getSize(), DataVector(transformation.getNrows()));
  for (size_t c = 0; c < input.getSize(); c++) {
    newDist[c] = Tools::mult(transformation, input.getVectors()[c]);
  }
  return FixedDistribution(input.getSize(), transformation.getNrows(), newDist);
}

sgpp::base::DataMatrix centerMean(sgpp::base::VectorDistribution& input) {
  std::vector<double> means(input.getDimensions());
  for (size_t d = 0; d < input.getDimensions(); ++d) {
    double mean = 0;
    for (size_t c = 0; c < input.getSize(); c++) {
      mean += input.getVectors()[c][d];
    }
    mean /= static_cast<double>(input.getSize());
    means[d] = mean;
  }

  sgpp::base::DataMatrix newVectors = input.getAsDataMatrix();
  for (size_t c = 0; c < input.getSize(); c++) {
    for (size_t d = 0; d < input.getDimensions(); ++d) {
      newVectors.set(c, d, newVectors.get(c, d) - means[d]);
    }
  }
  return std::move(newVectors);
}


sgpp::base::PcaInfo sgpp::base::PcaCovarianceSolver::solve(DataMatrix& matrix) {
  Eigen::MatrixXd b = Tools::toEigen(matrix);
  Eigen::MatrixXd c = (b.transpose() * b) * (1.0 / static_cast<double>(matrix.getNrows() - 1));
  PcaInfo i;
  i.basis = sgpp::base::DataMatrix(matrix.getNcols(), matrix.getNcols());
  i.eigenValues = sgpp::base::DataVector(matrix.getNcols());
  Tools::svd(c, i.basis, i.eigenValues);
  return i;
}


sgpp::base::PcaIterativeSolver::PcaIterativeSolver(size_t iterations, uint64_t seed) : iterations(iterations), seed(seed) {
}

sgpp::base::PcaInfo sgpp::base::PcaIterativeSolver::solve(DataMatrix& matrix) {
  size_t dimension = matrix.getNcols();
  sgpp::base::DataMatrix random(dimension, dimension);
  RandomOrientationDistribution orientations(dimension, seed, dimension, 1);
  for (size_t d = 0; d < dimension; ++d) {
    random.setColumn(d, orientations.getVectors()[d]);
  }

  Eigen::MatrixXd x = Tools::toEigen(matrix);
  Eigen::MatrixXd r = Tools::toEigen(random);
  Eigen::MatrixXd s(dimension, dimension);
  Eigen::MatrixXd e(dimension, dimension);
  for (size_t it = 0; it < iterations; ++it) {
    s = x.transpose() * (x * r);
    e = r.transpose() * s;
    r = s;
    for (size_t c = 0; c < dimension; c++) {
      r.col(c).normalize();
    }
  }

  
  PcaInfo i;
  i.basis = Tools::fromEigen(r);
  i.eigenValues = sgpp::base::DataVector(matrix.getNcols());
  for (size_t c = 0; c < dimension; c++) {
    i.eigenValues.set(c, e(c, c));
  }

  return i;
}

sgpp::base::PcaReducer::PcaReducer(std::shared_ptr<PcaSolver> solver) : solver(solver) {
}

const double MIN_EIGEN_VALUE = std::pow(10.0, -5.0);

sgpp::base::PcaInfo sgpp::base::PcaReducer::evaluate(VectorDistribution& input) {
  DataMatrix dist = centerMean(input);
  PcaInfo i = solver->solve(dist);

  size_t dimension = input.getDimensions();
  double sum = 0;
  i.activeComponentsCount = dimension;
  for (size_t d = 0; d < i.eigenValues.size(); ++d) {
    if (std::abs(i.eigenValues[d]) < MIN_EIGEN_VALUE) {
      i.activeComponentsCount = d;
      break;
      }
    sum += i.eigenValues[d];
  }

  i.principalAxes = i.basis;
  i.principalAxes.resizeRowsCols(dimension, i.activeComponentsCount);
  i.eigenValues.resize(i.activeComponentsCount);

  i.varianceShares = DataVector(i.activeComponentsCount);
  i.singularValues = DataVector(i.activeComponentsCount);
  i.loadings = DataMatrix(dimension, i.activeComponentsCount);
  for (size_t d = 0; d < i.activeComponentsCount; ++d) {
    i.varianceShares[d] = i.eigenValues[d] / static_cast<double>(sum);
    i.singularValues[d] = std::sqrt(i.eigenValues[d] * static_cast<double>(input.getSize() - 1));
    DataVector v(dimension);
    i.principalAxes.getColumn(d, v);
    v.mult(std::sqrt(i.eigenValues[d]));
    i.loadings.setColumn(d, v);
  }
  return i;
}