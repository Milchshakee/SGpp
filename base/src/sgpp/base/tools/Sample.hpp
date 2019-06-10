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
#include "sgpp/optimization/function/vector/VectorFunction.hpp"
#include "sgpp/datadriven/operation/hash/OperationPiecewiseConstantRegression/OperationPiecewiseConstantRegression.hpp"

namespace sgpp {
namespace base {

template <class K, class T>
class Sample {
 public:
  Sample() = default;

  Sample(const std::vector<K>& keys, std::function<T(K&)>& func) : keys(keys), values(keys.size()) {
    for (size_t i = 0; i < keys.size(); i++) {
      values[i] = func(keys[i]);
    }
  }

  Sample(const std::vector<K>& vectors, const std::vector<T>& values)
      : keys(vectors), values(values) {}

  const T& getValue(const K& key) {
    for (size_t i = 0; i < getSize(); i++) {
      if (key == keys[i]) {
        return values[i];
      }
    }
    throw std::invalid_argument();
  }
  size_t getSize() const { return values.size(); }
  size_t getDimensions() const { return keys.empty() ? 0 : keys[0].getSize(); }
  const std::vector<K>& getKeys() const { return keys; }
  const std::vector<T>& getValues() const { return values; }

 protected:
  std::vector<K> keys;
  std::vector<T> values;
};

template <class T>
using PointSample = Sample<DataVector, T>;

template <class T>
class GridSample : public PointSample<T> {
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
      : keys(grid.getSize()), values(values) {
    GridDistribution d(*grid);
    keys = d.getVectors();
  }

  const Grid& getGrid() const { return *grid; }

 private:
  std::shared_ptr<Grid> grid;
};

namespace SampleHelper {
template <class T>
PointSample<T> sampleDistribution(VectorDistribution& dist, std::function<T(DataVector&)>& func) {
  std::vector<double> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    values[i] = func(dist.getVectors()[i]);
  }
  return PointSample<T>(dist.getVectors(), values);
}

PointSample<double> sampleScalarFunction(VectorDistribution& dist,
                                         optimization::ScalarFunction& func) {
  std::vector<double> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    values[i] = func.eval(dist.getVectors()[i]);
  }
  return PointSample<double>(dist.getVectors(), values);
}

GridSample<double> sampleGrid(std::shared_ptr<Grid>& grid, optimization::ScalarFunction& func) {
  GridDistribution d(*grid);
  PointSample<double> s = sampleScalarFunction(d, func);
  return std::move(GridSample<double>(grid, s.getValues()));
}

PointSample<DataVector> sampleVectorFunction(VectorDistribution& dist,
                                             optimization::VectorFunction& func) {
  std::vector<DataVector> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    func.eval(dist.getVectors()[i], values[i]);
  }
  return PointSample<DataVector>(dist.getVectors(), values);
}

GridSample<DataVector> sampleGrid(std::shared_ptr<Grid>& grid,
                                   optimization::VectorFunction& func) {
  GridDistribution d(*grid);
  PointSample<DataVector> s = sampleVectorFunction(d, func);
  return std::move(GridSample<DataVector>(grid, s.getValues()));
}

  
  DataVector doHierarchisation(GridSample<double> sample) {
  DataVector v(sample.getValues());
  std::unique_ptr<sgpp::base::OperationHierarchisation>(
      sgpp::op_factory::createOperationHierarchisation(const_cast<Grid&>(sample.getGrid())))
      ->doHierarchisation(v);
  return std::move(v);
  }

}  // namespace SampleHelper

}  // namespace base
}  // namespace sgpp

#endif
