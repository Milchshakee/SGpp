#include "sgpp/base/tools/VectorDistribution.hpp"
#include "sgpp/base/tools/dimension/PcaReducer.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
  sgpp::base::DataMatrix data = sgpp::base::DataMatrix::fromString("[[0.3, 0.25],[-0.75, -0.5375],[0.9375, 0.9375],[-0.4375, -0.75],[0.75, 0.4375],[-1.0, -0.9375],[0.4375, 0.5625],[-0.5625, -0.5625],[0.5625, 0.65]]");
  //sgpp::base::DataMatrix data = sgpp::base::DataMatrix::fromString("[[-1, -1],[0.5, 0.5]]");
  auto dist = sgpp::base::FixedDistribution(data);
  auto cutter = sgpp::base::PcaVarianceCutter(0.95);
  auto reducer = sgpp::base::PcaReducer(10000, std::mt19937_64::default_seed);

  auto info = reducer.evaluate(dist);
      for (size_t d = 0; d < dist.getDimensions(); ++d) {
    std::cout << "dimension: " + std::to_string(d) << std::endl;
        sgpp::base::DataVector col(2);
    info.eigenVectors.getColumn(d, col);
    std::cout << "eigen vector: " + col.toString() << std::endl;
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
