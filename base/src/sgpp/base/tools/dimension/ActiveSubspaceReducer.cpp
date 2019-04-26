#include "ActiveSubspaceReducer.hpp"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

ActiveSubspaceReducer::MatrixComputationStrategy::MatrixComputationStrategy(
    size_t dimensions, GradientGenerationStrategy& s)
    : dimensions(dimensions), gradient(s) {}

ActiveSubspaceReducer::MCStrategy::MCStrategy(size_t dimensions, GradientGenerationStrategy& s,
                                              VectorDistribution& distribution, size_t samples)
    : MatrixComputationStrategy(dimensions, s), distribution(distribution), samples(samples) {}

std::unique_ptr<sgpp::base::DataMatrix> ActiveSubspaceReducer::MCStrategy::compute() {
  auto result = std::make_unique<sgpp::base::DataMatrix>(dimensions, dimensions);

  for (size_t i = 0; i < samples; ++i) {
    sgpp::base::DataVector sample = distribution();
    sgpp::base::DataVector sampleGradient = gradient.gradientAt(sample);
    sgpp::base::DataMatrix temp(dimensions, dimensions);
    for (size_t d = 0; d < dimensions; ++d) {
      sgpp::base::DataVector col = sampleGradient;
      col.mult(sampleGradient[d]);
      temp.setColumn(d, col);
    }
    result->add(temp);
  }

  result->mult(1.0 / samples);
  return std::move(result);
}

std::unique_ptr<sgpp::optimization::VectorFunction> ActiveSubspaceReducer::reduce(
    sgpp::optimization::VectorFunction& input)
{
  std::unique_ptr<sgpp::base::DataMatrix> c = matrix.compute();

  size_t n = lhsMatrix.getNrows();

  gsl_matrix_view m = gsl_matrix_view_array(lhsMatrix.getPointer(), n,
                                            n);  // Create GSL matrix view for decomposition

  auto q = std::unique_ptr<gsl_matrix>{gsl_matrix_alloc(n, n)};  // Stores the eigenvectors
  auto e = std::unique_ptr<gsl_vector>{gsl_vector_alloc(n)};     // Stores the eigenvalues

  gsl_eigen_symmv_workspace* ws = gsl_eigen_symmv_alloc(n);
  gsl_eigen_symmv(&m.matrix, e.get(), q.get(), ws);
  gsl_eigen_symmv_free(ws);

  // Create an (n+1)*n matrix to store eigenvalues and -vectors:
  lhsMatrix = sgpp::base::DataMatrix(n + 1, n);

  for (size_t r = 0; r < n; r++) {
    for (size_t c = 0; c < n; c++) {
      lhsMatrix.set(r, c, gsl_matrix_get(q.get(), r, c));
    }
  }
  for (size_t c = 0; c < n; c++) {
    lhsMatrix.set(n, c, gsl_vector_get(e.get(), c));
  }

  sgpp::base::DataMatrix eigenVectorMatrix; //TODO
  sgpp::base::DataMatrix eigenValueMatrix; //TODO

}
