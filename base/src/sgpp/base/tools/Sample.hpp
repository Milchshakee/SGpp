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

template <class T, class K = DataVector>
class Sample {
 public:
  Sample() = default;

  Sample(const std::vector<K>& keys, std::function<T(K&)>& func)
      : keys(keys), values(keys.size())
  {
    for (size_t i = 0; i < keys.size(); i++) {
      values[i] = func(keys[i]);
    }
  }

  Sample(const std::vector<K>& vectors, const std::vector<T>& values)
      : keys(vectors), values(values) {}

  size_t getSize() const { return values.size(); }
  size_t getDimensions() const { return keys.empty() ? 0 : keys[0].getSize(); }
  const std::vector<K>& getKeys() const { return keys; }
  const std::vector<T>& getValues() const { return values; }

 protected:
  std::vector<K> keys;
  std::vector<T> values;
};

  template <class T>
class GridSample : public Sample<T> {
 public:
  GridSample() = default;

  GridSample(std::shared_ptr<Grid>& grid, std::function<T(DataVector&)>& func)
      : keys(grid->getSize()), values(grid->getSize()) {
    GridDistribution d(*grid);
    for (size_t i = 0; i < d.getSize(); i++) {
      values[i] = func(d.getVectors()[i]);
    }
  }

  GridSample(std::shared_ptr<Grid>& grid, const std::vector<T>& values)
      : keys(grid.getSize()), values(values)
  {
    GridDistribution d(*grid);
    keys = d.getVectors();
  }

    const Grid& getGrid() const { return *grid; }

 private:
  std::shared_ptr<Grid> grid;
  };

namespace Sampler {

template <class T>
Sample<T> sampleDistribution(VectorDistribution& dist, std::function<T(DataVector&)>& func) {
  std::vector<double> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    values[i] = func(dist.getVectors()[i]);
  }
  return Sample<T>(dist.getVectors(), values);
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
