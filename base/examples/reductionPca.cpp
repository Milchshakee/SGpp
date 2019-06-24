#include <iostream>
#include <sgpp/base/tools/VectorDistribution.hpp>
#include <sgpp/base/tools/dimension/PcaReducer.hpp>

int main(int argc, char* argv[]) {
  // Create the data set.
  sgpp::base::DataMatrix data = sgpp::base::DataMatrix::fromString(
      "[[0, 0, 0],[1, 1, 2],[1,0,1],[0,1,1],[0.5,0.5,1]]");
  // Create a wrapper of the raw data for easier usage
  auto dist = sgpp::base::FixedDistribution(data);

  // Create the solver used to calculate the principal axes and their eigen values
  auto solver = std::make_shared<sgpp::base::PcaCovarianceSolver>();
  // auto solver = std::make_shared<sgpp::base::PcaIterativeSolver>(1000, std::mt19937_64::default_seed);

  // Create the reducer.
  auto reducer = sgpp::base::PcaReducer(solver);

  // Use the reducer to first evaluate the data
  sgpp::base::PcaInfo info = reducer.evaluate(dist);

  // Print out all the information that the reducer has gathered from the data

  //First, loop over all active principal components (principal directions with an eigen value bigger than 0)
  for (size_t d = 0; d < info.activeComponentsCount; ++d) {
    std::cout << "dimension: " + std::to_string(d) << std::endl;
    sgpp::base::DataVector pa(3);
    info.principalAxes.getColumn(d, pa);
    std::cout << "principal direction: " + pa.toString() << std::endl;
    sgpp::base::DataVector l(3);
    info.loadings.getColumn(d, l);
    std::cout << "loading: " + l.toString() << std::endl;
    std::cout << "eigen value: " + std::to_string(info.eigenValues[d]) << std::endl;
    std::cout << "singular value: " + std::to_string(info.singularValues[d]) << std::endl;
    std::cout << "share of variance: " + std::to_string(info.varianceShares[d]) << std::endl;
    std::cout << std::endl;
  }

  // If info.activeComponentsCount != dimensions, we know that there are some inactive principal directions with an eigen value of 0
  for (size_t d = info.activeComponentsCount; d < dist.getDimensions(); ++d) {
    std::cout << "dimension: " + std::to_string(d) << " (inactive)" << std::endl;
    sgpp::base::DataVector pa(3);
    info.basis.getColumn(d, pa);
    std::cout << "principal direction: " + pa.toString() << std::endl;
    std::cout << "eigen value: 0" << std::endl;
    std::cout << std::endl;
  }

  // Create a cutter to reduce the dimensions
  // In this case, we want only want to preserve the most dominant principal axes that comprise at
  // least 95% of the data variance.
  auto cutter = sgpp::base::PcaVarianceCutter(0.95);
  // Alternatively, we can also reduce the dimensions to a fixed parameter without caring about the
  // covered variance of that reduction.
  // auto cutter = sgpp::base::PcaFixedCutter(1);

  // Use the cutter to remove unwanted dimensions
  sgpp::base::PcaResult result = cutter.cut(dist, info);

  // Print information of the result object.
  std::cout << "transformation matrix used to reduce the dataset:\n" +
                   result.transformation.toString()
            << std::endl;
  std::cout << "actual covered variance after reduction: " + std::to_string(result.coveredVariance)
            << std::endl;
  std::cout << std::endl;

  // Apply the result to the actual data set.
  // In this step, the data vectors are transformed using the transformation matrix from the result.
  auto reduced = result.apply(dist);
  sgpp::base::DataMatrix reducedData = reduced.getAsDataMatrix();

  // Print information about the reduced data set.
  std::cout << "dimensions after reduction: " + std::to_string(reduced.getDimensions()) << std::endl;
  std::cout << "reduced data:\n" + reducedData.toString() << std::endl;

  return 0;
}
