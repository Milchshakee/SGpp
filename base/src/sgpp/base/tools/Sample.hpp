// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <vector>
#include "VectorDistribution.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/optimization/function/scalar/ScalarFunction.hpp"

namespace sgpp {
namespace base {

template <class T>
class Sample {
 public:
  Sample() = default;

  Sample(VectorDistribution& dist, std::function<T(DataVector&)>& func)
      : vectors(dist.getVectors()), values(dist.getSize()) {
    for (size_t i = 0; i < dist.getSize(); i++) {
      values[i] = func(dist.getVectors()[i]);
    }
  }

  Sample(const std::vector<DataVector>& vectors, const std::vector<T>& values)
      : vectors(vectors), values(values) {}

  size_t getSize() const { return values.size(); }
  size_t getDimensions() const { return vectors.empty() ? 0 : vectors[0].getSize(); }
  const std::vector<DataVector>& getVectors() const { return vectors; }
  const std::vector<T>& getValues() const { return values; }

 private:
  std::vector<DataVector> vectors;
  std::vector<T> values;
};

namespace Sampler {
template <class T>
Sample<T> sampleGrid(Grid& grid, std::function<T(DataVector)>& func) {
  GridDistribution d(grid);
  return Sample<T>(d, func);
}

Sample<double> sampleScalarFunction(VectorDistribution& dist,
                                           optimization::ScalarFunction& func) {
  std::vector<double> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    values[i] = func.eval(dist.getVectors()[i]);
  }
  return Sample<double>(dist.getVectors(), values);
}

Sample<double> sampleGrid(Grid& grid, optimization::ScalarFunction& func) {
  GridDistribution d(grid);
  return sampleScalarFunction(d, func);
}

Sample<DataVector> sampleVectorFunction(VectorDistribution& dist,
                                               optimization::VectorFunction& func) {
  std::vector<DataVector> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    func.eval(dist.getVectors()[i], values[i]);
  }
  return Sample<DataVector>(dist.getVectors(), values);
}

Sample<DataVector> sampleGrid(Grid& grid, optimization::VectorFunction& func) {
  GridDistribution d(grid);
  return sampleVectorFunction(d, func);
}

}  // namespace Sampler

}  // namespace base
}  // namespace sgpp

#endif
