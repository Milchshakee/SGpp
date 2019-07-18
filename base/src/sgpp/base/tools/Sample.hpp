// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include <vector>
#include "OperationQuadratureMC.hpp"
#include "VectorDistribution.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/optimization/function/scalar/ScalarFunction.hpp"
#include "sgpp/optimization/function/vector/VectorFunction.hpp"

namespace sgpp {
namespace base {

template <class K, class T>
class Sample {
 public:
  Sample() = default;

  Sample(const std::vector<K>& keys, std::function<T(const K&)>& func)
      : keys(keys), values(keys.size()) {
    for (size_t i = 0; i < keys.size(); i++) {
      values[i] = func(keys[i]);
    }
  }

  Sample(const std::map<K, T> map) : keys(map.size()), values(map.size()) {
    size_t counter = 0;
    for (auto i = map.begin(); i != map.end(); ++i) {
      keys[counter] = i->first;
      values[counter] = i->second;
      counter++;
    }
  }

  Sample(const std::vector<K>& vectors, const std::vector<T>& values)
      : keys(vectors), values(values) {}

  const T& getValue(const K& key) const {
    for (size_t i = 0; i < getSize(); i++) {
      if (key == keys[i]) {
        return values[i];
      }
    }
    throw std::invalid_argument("Value not found");
  }
  size_t getSize() const { return values.size(); }
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
      : Sample<DataVector, T>(keys, func) {}

  PointSample(const std::vector<DataVector>& keys, const std::vector<T>& values)
      : Sample<DataVector, T>(keys, values) {}

  size_t getDimensions() const {
    return Sample<DataVector, T>::keys.empty() ? 0 : Sample<DataVector, T>::keys[0].getSize();
  }
};

template <class T>
class GridSample : public PointSample<T> {
 public:
  GridSample() = default;

  GridSample(std::shared_ptr<Grid>& grid, std::function<T(const DataVector&)>& func) : grid(grid) {
    PointSample<T>::keys = std::vector<DataVector>(grid->getSize());
    PointSample<T>::values = std::vector<T>(grid->getSize());
    GridDistribution d(*grid);
    for (size_t i = 0; i < d.getSize(); i++) {
      PointSample<T>::keys[i] = d.getVectors()[i];
      PointSample<T>::values[i] = func(d.getVectors()[i]);
    }
  }

  GridSample(std::shared_ptr<Grid>& grid, const std::vector<T>& values) : grid(grid) {
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
      : GridSample<double>(grid, func), valuesView(values.data(), values.size()) {}

  SGridSample(std::shared_ptr<Grid>& grid, optimization::ScalarFunction& func) {
    GridSample<double>::grid = grid;
    keys = std::vector<DataVector>(grid->getSize());
    values = std::vector<double>(grid->getSize());

    GridDistribution d(*grid);
    for (size_t i = 0; i < d.getSize(); i++) {
      keys[i] = d.getVectors()[i];
      values[i] = func.eval(d.getVectors()[i]);
    }
    valuesView = DataVector(values.data(), values.size());
  }

  SGridSample(std::shared_ptr<Grid>& grid, const std::vector<double>& values)
      : GridSample<double>(grid, values),
        valuesView(GridSample<double>::values.data(), GridSample<double>::values.size()) {}

  void hierarchise() {
    std::unique_ptr<sgpp::base::OperationHierarchisation>(
        sgpp::op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(valuesView);
  }

  double eval(const DataVector& point) {
    std::unique_ptr<sgpp::base::OperationEval> op(sgpp::op_factory::createOperationEval(*grid));
    return op->eval(valuesView, point);
  }

  double quadrature() {
    std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
        sgpp::op_factory::createOperationQuadrature(*grid));
    double res = opQ->doQuadrature(valuesView);
    return res;
  }

  double mcQuadrature(size_t paths) {
    sgpp::base::OperationQuadratureMC opMC(*grid, paths);
    return opMC.doQuadrature(valuesView);
  }

  double mcL2Error(FUNC f, void* clientdata, size_t paths) {
    sgpp::base::OperationQuadratureMC opMC(*grid, paths);
    return opMC.doQuadratureL2Error(f, clientdata, valuesView);
  }

  double mcL2Error(optimization::ScalarFunction& f, size_t paths) {
    sgpp::base::OperationQuadratureMC opMC(*grid, paths);
    return opMC.doQuadratureL2Error(f, valuesView);
  }

  const DataVector& getValuesView() const { return valuesView; }

 private:
  DataVector valuesView;
};

namespace SampleHelper {
template <class T>
PointSample<T> sampleDistribution(VectorDistribution& dist,
                                  std::function<T(const DataVector&)>& func) {
  std::vector<double> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    values[i] = func(dist.getVectors()[i]);
  }
  return PointSample<T>(dist.getVectors(), values);
}

inline PointSample<double> sampleScalarFunction(VectorDistribution& dist,
                                                optimization::ScalarFunction& func) {
  std::vector<double> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    values[i] = func.eval(dist.getVectors()[i]);
  }
  return PointSample<double>(dist.getVectors(), values);
}

inline PointSample<DataVector> sampleVectorFunction(VectorDistribution& dist,
                                                    optimization::VectorFunction& func) {
  std::vector<DataVector> values(dist.getSize());
  for (size_t i = 0; i < dist.getSize(); i++) {
    func.eval(dist.getVectors()[i], values[i]);
  }
  return PointSample<DataVector>(dist.getVectors(), values);
}

inline GridSample<DataVector> sampleGrid(std::shared_ptr<Grid>& grid,
                                         optimization::VectorFunction& func) {
  GridDistribution d(*grid);
  PointSample<DataVector> s = sampleVectorFunction(d, func);
  return std::move(GridSample<DataVector>(grid, s.getValues()));
}

}  // namespace SampleHelper
}  // namespace base
}  // namespace sgpp

#endif
