// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/exception/tool_exception.hpp>
#include <sgpp/base/tools/DistributionSample.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

namespace sgpp {
namespace base {

template <class K, class T, class V = std::vector<T>>
class Sample {
public:
  Sample() = default;

  Sample(const std::vector<K>& keys, std::function<T(const K&)>& func)
    : keys(keys),
      values(keys.size()) {
    if (keys.size() != values.size()) {
      throw tool_exception("Key size and value size do not match");
      }
    for (size_t i = 0; i < keys.size(); i++) {
      values[i] = func(keys[i]);
    }
  }

  Sample(const std::map<K, T> map)
    : keys(map.size()), values(map.size()) {
    if (keys.size() != values.size()) {
      throw tool_exception("Key size and value size do not match");
    }
    size_t counter = 0;
    for (auto i = map.begin(); i != map.end(); ++i) {
      keys[counter] = i->first;
      values[counter] = i->second;
      counter++;
    }
  }

  Sample(const std::vector<K>& vectors, const std::vector<T>& values)
    : keys(vectors), values(values) {
    if (keys.size() != values.size()) {
      throw tool_exception("Key size and value size do not match");
    }
  }

  void erase(size_t index)
  {
    keys.erase(keys.begin() + index);
    values.erase(values.begin() + index);
  }

  T& getValue(const K& key) {
    for (size_t i = 0; i < getSize(); i++) {
      if (key == keys[i]) {
        return values[i];
      }
    }
    throw std::invalid_argument("Value not found");
  }

  const T& getValue(const K& key) const {
    for (size_t i = 0; i < getSize(); i++) {
      if (key == keys[i]) {
        return values[i];
      }
    }
    throw std::invalid_argument("Value not found");
  }

  size_t getSize() const { return values.size(); }
  std::vector<K>& getKeys() { return keys; }
  V& getValues() { return values; }
  const std::vector<K>& getKeys() const { return keys; }
  const V& getValues() const { return values; }

protected:
  std::vector<K> keys;
  V values;
};

template <class T, class V = std::vector<T>>
class PointSample : public Sample<DataVector, T, V> {
public:
  PointSample() = default;

  PointSample(const std::vector<DataVector>& keys, std::function<T(const DataVector&)>& func)
    : Sample<DataVector, T, V>(keys, func) {
  }

  PointSample(const std::vector<DataVector>& keys, const std::vector<T>& values)
    : Sample<DataVector, T, V>(keys, values) {
  }

  size_t getDimensions() const {
    return Sample<DataVector, T, V>::keys.empty() ? 0 : Sample<DataVector, T, V>::keys[0].getSize();
  }
};

template <class T, class V = std::vector<T>>
class GridSample : public PointSample<T,V> {
public:
  GridSample() = default;

  GridSample(std::shared_ptr<Grid>& grid, std::function<T(const DataVector&)>& func)
    : grid(grid) {
    PointSample<T, V>::keys = std::vector<DataVector>(grid->getSize());
    PointSample<T, V>::values = V(grid->getSize());
    DistributionSample d(*grid);
    for (size_t i = 0; i < d.getSize(); i++) {
      PointSample<T,V>::keys[i] = d.getVectors()[i];
      PointSample<T,V>::values[i] = func(d.getVectors()[i]);
    }
  }

  GridSample(std::shared_ptr<Grid>& grid, const V& values) : grid(grid) {
    DistributionSample d(*grid);
    PointSample<T, V>::keys = d.getVectors();
    PointSample<T, V>::values = values;
  }

  const Grid& getGrid() const { return *grid; }

protected:
  std::shared_ptr<Grid> grid;
};

class SGridSample : public GridSample<double, DataVector> {
public:
  SGridSample() = default;

  SGridSample(std::shared_ptr<Grid>& grid, std::function<double(const DataVector&)>& func)
    : GridSample<double,DataVector>(grid, func), hierarchised(false) {
  }

  SGridSample(std::shared_ptr<Grid>& grid, ScalarFunction& func) : hierarchised(false) {
    GridSample<double, DataVector>::grid = grid;
    keys = std::vector<DataVector>(grid->getSize());
    values = DataVector(values.data(), values.size());

    DistributionSample d(*grid);
    for (size_t i = 0; i < d.getSize(); i++) {
      keys[i] = d.getVectors()[i];
      values[i] = func.eval(d.getVectors()[i]);
    }
  }

  SGridSample(std::shared_ptr<Grid>& grid, const DataVector& values)
    : GridSample<double,DataVector>(grid, values),
        hierarchised(false) {
  }

  void hierarchise() {
    if (hierarchised) {
      throw tool_exception("Data is already hierarchised");
      }

    std::unique_ptr<OperationHierarchisation>(
          op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(values);
    hierarchised = true;
  }

    void dehierarchise() {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    std::unique_ptr<OperationHierarchisation>(op_factory::createOperationHierarchisation(*grid))
        ->doDehierarchisation(values);
    hierarchised = false;
  }

  double eval(const DataVector& point) const {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    std::unique_ptr<OperationEval> op(op_factory::createOperationEval(*grid));
    return op->eval(values, point);
  }

  double quadrature() const {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    std::unique_ptr<OperationQuadrature> opQ(
      op_factory::createOperationQuadrature(*grid));
    double res = opQ->doQuadrature(const_cast<DataVector&>(values));
    return res;
  }

  double mcQuadrature(size_t paths) {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    OperationQuadratureMC opMC(*grid, paths);
    return opMC.doQuadrature(values);
  }

  void setHierarchised(bool h) { hierarchised = h; }

  bool isHierarchised() const { return hierarchised; }

private:
  bool hierarchised;
};

namespace SampleHelper {

inline PointSample<double> sampleScalarFunction(DistributionSample& dist, ScalarFunction& func) {
  std::vector<double> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    values[i] = func.eval(dist.getVectors()[i]);
  }
  return PointSample<double>(dist.getVectors(), values);
}

inline PointSample<DataVector> sampleVectorFunction(DistributionSample& dist,
                                                    VectorFunction& func) {
  std::vector<DataVector> values(dist.getSize(), DataVector(dist.getDimensions()));
  for (size_t i = 0; i < dist.getSize(); i++) {
    func.eval(dist.getVectors()[i], values[i]);
  }
  return PointSample<DataVector>(dist.getVectors(), values);
}

inline GridSample<DataVector> sampleGrid(std::shared_ptr<Grid>& grid, VectorFunction& func) {
  DistributionSample d(*grid);
  PointSample<DataVector> s = sampleVectorFunction(d, func);
  return std::move(GridSample<DataVector>(grid, s.getValues()));
}

  inline datadriven::Dataset fromPointSample(PointSample<double>& sample) {
  datadriven::Dataset d(sample.getSize(), sample.getDimensions());
  for (size_t i = 0; i < sample.getSize(); i++) {
    d.getData().setRow(i, sample.getKeys()[i]);
    d.getTargets()[i] = sample.getValues()[i];
  }
  return d;
  }
} // namespace SampleHelper
} // namespace base
} // namespace sgpp
