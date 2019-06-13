#include "sgpp/base/tools/VectorDistribution.hpp"
#include "sgpp/base/tools/dimension/PcaReducer.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
  sgpp::base::DataMatrix data = sgpp::base::DataMatrix::fromFile("pcaReduction.input.txt");
  auto dist = sgpp::base::FixedDistribution(data);
  auto cutter = sgpp::base::PcaVarianceCutter(0.95);
  auto reducer = sgpp::base::PcaReducer(1000, std::mt19937_64::default_seed);

  auto info = reducer.evaluate(dist);
      for (size_t d = 0; d < dist.getDimensions(); ++d) {
    std::cout << "dimension: " + d << std::endl;
    std::cout << "eigen vector: " + std::to_string(info.eigenVectors[d]) << std::endl;
    std::cout << "eigen value: " + std::to_string(info.eigenValues[d]) << std::endl;
    std::cout << "share of variance: " + std::to_string(info.varianceShares[d]) << std::endl;
  }


  auto result = cutter.cut(dist, info);
  auto reduced = result.apply(dist);
  sgpp::base::DataMatrix reducedData = reduced.getAsDataMatrix();
  std::cout << "transformation matrix used to reduce the dataset: " + result.transformation.toString() << std::endl;
  std::cout << "actual covered variance after reduction: " + std::to_string(result.coveredVariance) << std::endl;
  std::cout << "dimensions after reduction: " + reduced.getDimensions() << std::endl;
  std::cout << "reduced data: " + reducedData.toString() << std::endl;

  return 0;
}
