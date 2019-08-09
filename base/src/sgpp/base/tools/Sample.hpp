// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>
#include <sgpp/base/tools/dist/GridDistribution.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/exception/tool_exception.hpp>
#include <sgpp/base/tools/OperationL2.hpp>

namespace sgpp {
namespace base {

template <class K, class T>
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
  std::vector<T>& getValues() { return values; }
  const std::vector<K>& getKeys() const { return keys; }
  const std::vector<T>& getValues() const { return values; }

protected:
  std::vector<K> keys;
  std::vector<T> values;
};

template <class T>
class PointSample : public Sample<DataVector, T> {
public:
  PointSample() = default;

  PointSample(const std::vector<DataVector>& keys, std::function<T(const DataVector&)>& func)
    : Sample<DataVector, T>(keys, func) {
  }

  PointSample(const std::vector<DataVector>& keys, const std::vector<T>& values)
    : Sample<DataVector, T>(keys, values) {
  }

  size_t getDimensions() const {
    return Sample<DataVector, T>::keys.empty() ? 0 : Sample<DataVector, T>::keys[0].getSize();
  }
};

template <class T>
class GridSample : public PointSample<T> {
public:
  GridSample() = default;

  GridSample(std::shared_ptr<Grid>& grid, std::function<T(const DataVector&)>& func)
    : grid(grid) {
    PointSample<T>::keys = std::vector<DataVector>(grid->getSize());
    PointSample<T>::values = std::vector<T>(grid->getSize());
    GridDistribution d(*grid);
    for (size_t i = 0; i < d.getSize(); i++) {
      PointSample<T>::keys[i] = d.getVectors()[i];
      PointSample<T>::values[i] = func(d.getVectors()[i]);
    }
  }

  GridSample(std::shared_ptr<Grid>& grid, const std::vector<T>& values)
    : grid(grid) {
    PointSample<T>::keys = std::vector<DataVector>(grid->getSize());
    PointSample<T>::values = std::vector<T>(values);
    GridDistribution d(*grid);
    PointSample<T>::keys = d.getVectors();
  }

  const Grid& getGrid() const { return *grid; }

protected:
  std::shared_ptr<Grid> grid;
};

class SGridSample : public GridSample<double> {
public:
  SGridSample() = default;

  SGridSample(std::shared_ptr<Grid>& grid, std::function<double(const DataVector&)>& func)
    : GridSample<double>(grid, func),
      valuesDataVector(values.data(), values.size()), hierarchised(false) {
  }

  SGridSample(std::shared_ptr<Grid>& grid, ScalarFunction& func) : hierarchised(false) {
    GridSample<double>::grid = grid;
    keys = std::vector<DataVector>(grid->getSize());
    values = std::vector<double>(grid->getSize());

    GridDistribution d(*grid);
    for (size_t i = 0; i < d.getSize(); i++) {
      keys[i] = d.getVectors()[i];
      values[i] = func.eval(d.getVectors()[i]);
    }
    valuesDataVector = DataVector(values.data(), values.size());
  }

  SGridSample(std::shared_ptr<Grid>& grid, const std::vector<double>& values)
    : GridSample<double>(grid, values),
        valuesDataVector(GridSample<double>::values.data(), GridSample<double>::values.size()),
        hierarchised(false) {
  }

  void hierarchise() {
    if (hierarchised) {
      throw tool_exception("Data is already hierarchised");
      }

    std::unique_ptr<OperationHierarchisation>(
          op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(valuesDataVector);
    sync();
    hierarchised = true;
  }

    void dehierarchise() {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    std::unique_ptr<OperationHierarchisation>(op_factory::createOperationHierarchisation(*grid))
        ->doDehierarchisation(valuesDataVector);
    sync();
    hierarchised = false;
  }

  double eval(const DataVector& point) const {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    std::unique_ptr<OperationEval> op(op_factory::createOperationEval(*grid));
    return op->eval(valuesDataVector, point);
  }

  double quadrature() const {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    std::unique_ptr<OperationQuadrature> opQ(
      op_factory::createOperationQuadrature(*grid));
    double res = opQ->doQuadrature(const_cast<DataVector&>(valuesDataVector));
    return res;
  }

  double mcQuadrature(size_t paths) {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    OperationQuadratureMC opMC(*grid, paths);
    return opMC.doQuadrature(valuesDataVector);
  }

  double mcL2Error(FUNC f, void* clientdata, size_t paths) {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    OperationQuadratureMC opMC(*grid, paths);
    return opMC.doQuadratureL2Error(f, clientdata, valuesDataVector);
  }

  double mcL2Error(ScalarFunction& f, size_t paths) {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

    OperationQuadratureMC opMC(*grid, paths);
    return opMC.doQuadratureL2Error(f, valuesDataVector);
  }

    double l2Norm() const {
    if (!hierarchised) {
      throw tool_exception("Data is not hierarchised");
    }

      DataVector v(values.size());
    double res = 0;
    for (size_t i = 0; i < getSize(); i++) {
      v[i] = values[i] * values[i];
    }
    OperationL2 op(grid->getStorage());
    return op.calculateL2Norm(valuesDataVector);
  }

  bool isHierarchised() const { return hierarchised; }

  const DataVector& getValuesDataVector() const { return valuesDataVector; }

  DataVector& getValuesDataVector() { return valuesDataVector; }

private:
  void sync() { values = std::vector<double>(valuesDataVector.begin(), valuesDataVector.end()); }

  DataVector valuesDataVector;
  bool hierarchised;
};

namespace SampleHelper {

inline PointSample<double> sampleScalarFunction(VectorDistribution& dist, ScalarFunction& func) {
  std::vector<double> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    values[i] = func.eval(dist.getVectors()[i]);
  }
  return PointSample<double>(dist.getVectors(), values);
}

inline PointSample<DataVector> sampleVectorFunction(VectorDistribution& dist,
                                                    VectorFunction& func) {
  std::vector<DataVector> values(dist.getSize(), DataVector(dist.getDimensions()));
  for (size_t i = 0; i < dist.getSize(); i++) {
    func.eval(dist.getVectors()[i], values[i]);
  }
  return PointSample<DataVector>(dist.getVectors(), values);
}

inline GridSample<DataVector> sampleGrid(std::shared_ptr<Grid>& grid, VectorFunction& func) {
  GridDistribution d(*grid);
  PointSample<DataVector> s = sampleVectorFunction(d, func);
  return std::move(GridSample<DataVector>(grid, s.getValues()));
}

} // namespace SampleHelper
} // namespace base
} // namespace sgpp
