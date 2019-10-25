

#include <sgpp/base/tools/dimension/AsMcReducer.hpp>

namespace sgpp {
namespace base {

VectorFunction& AsMcResult::getTransformationFunction() { return projection.getFunction(); }

AsMcFixedCutter::AsMcFixedCutter(size_t n, GridType type, level_t level, const DataVector& mean)
    : FixedCutter<sgpp::base::AsMcInput, sgpp::base::AsInfo, sgpp::base::AsMcResult>(n),
      AsMcCutter(type, level, mean) {}

AsMcErrorRuleCutter::AsMcErrorRuleCutter(ErrorRule& r, double maxError, GridType type,
                                         level_t level, const DataVector& mean)
    : ErrorRuleCutter<sgpp::base::AsMcInput, sgpp::base::AsInfo, sgpp::base::AsMcResult>(r,
                                                                                         maxError),
      AsMcCutter(type, level, mean) {}

AsMcResult AsMcErrorRuleCutter::cut(const AsMcInput& input, const AsInfo& info) {
  AsMcResult last(input, info.eigenVectors, input.function.getNumberOfParameters(), type, level, mean);
  for (size_t d = 1; d < input.function.getNumberOfParameters(); ++d) {
    AsMcResult result(input, info.eigenVectors, input.function.getNumberOfParameters() - d, type,
                      level, mean);
    double var = 1 - result.calculateRelativeError(r);
    if (var > minVariance) {
      last = result;
    } else {
      return last;
    }
  }
  return last;
}


ScalarFunction& AsMcResult::getOriginalFunction() { return *originalFunction; }

ScalarFunction& AsMcResult::getReducedFunctionSurrogate() { return evalFunc; }

SGridSample& AsMcResult::getReducedOutput() { return reduced; }


InputProjection& AsMcResult::getProjection() { return projection; }

PointSample<DataMatrix> AsMcReducer::fromGradientSample(const PointSample<DataVector>& gradients) {
  std::vector<DataMatrix> out(gradients.getSize(),
                              DataMatrix(gradients.getDimensions(), gradients.getDimensions()));
  for (size_t i = 0; i < gradients.getSize(); ++i) {
    sgpp::base::DataVector sampleGradient = gradients.getValues()[i];
    for (size_t d = 0; d < gradients.getDimensions(); ++d) {
      sgpp::base::DataVector col = sampleGradient;
      col.mult(sampleGradient[d]);
      out[i].setColumn(d, col);
    }
  }
  return PointSample<DataMatrix>(gradients.getKeys(), out);
}

PointSample<DataMatrix> AsMcReducer::fromFiniteDifferences(ScalarFunction& func,
                                                           VectorDistribution& v, double h) {
  std::vector<DataMatrix> samples(v.getSize(), DataMatrix(v.getDimensions(), v.getDimensions()));
  sgpp::base::DataVector sampleGradient(v.getDimensions());
  DataVector working;
  for (size_t i = 0; i < v.getSize(); ++i) {
    const DataVector& vec = v.getVectors()[i];
    working = vec;
    for (size_t d = 0; d < v.getDimensions(); ++d) {
      double hOffset = vec[d] + h;
      if (hOffset > 1.0) {
        hOffset -= 2 * h;
      }
      working[d] = hOffset;
      double val = (func.eval(vec) - func.eval(working)) / h;
      sampleGradient[d] = val;
      working[d] = vec[d];
    }

    for (size_t d = 0; d < v.getDimensions(); ++d) {
      sgpp::base::DataVector col = sampleGradient;
      col.mult(sampleGradient[d]);
      samples[i].setColumn(d, col);
    }
  }
  return PointSample<DataMatrix>(v.getVectors(), samples);
}

AsMcResult::AsMcResult(const AsMcInput& input, const DataMatrix& m, size_t n, GridType type,
                       level_t l, const DataVector& mean)
    : originalFunction(&input.function), projection(m, n, mean) {
  RegularGridConfiguration g;
  g.dim_ = n;
  g.level_ = l;
  g.type_ = type;
  std::shared_ptr<Grid> newGrid(Grid::createGrid(g));
  newGrid->getGenerator().regular(l);

  std::function<double(const DataVector&)> func = [this, &input](const DataVector& v) {
    DataVector newV;
    projection.inverse(v, newV);
    return input.function.eval(newV);
  };

  SGridSample newSample(newGrid, func);
  newSample.hierarchise();
  reduced = newSample;
  evalFunc = EvalFunction(reduced);
}

AsMcResult AsMcFixedCutter::cut(const AsMcInput& input, const AsInfo& info) {
  return AsMcResult(input, info.eigenVectors, n, type, level, mean);
}

AsMcIntervalCutter::AsMcIntervalCutter(size_t bootstrapSamples, GridType type, level_t level,
                                       const DataVector& mean)
    : AsMcCutter(type, level, mean), bootstrapSamples(bootstrapSamples) {}

AsInfo AsMcReducer::evaluate(AsMcInput& input) {
  size_t dimensions = input.samples.getDimensions();
  sgpp::base::DataMatrix matrix(dimensions, dimensions);
  for (size_t i = 0; i < input.samples.getSize(); ++i) {
    matrix.add(input.samples.getValues()[i]);
  }
  matrix.mult(1.0 / static_cast<double>(input.samples.getSize()));

  AsInfo i;
  i.eigenVectors = sgpp::base::DataMatrix(dimensions, dimensions);
  i.eigenValues = sgpp::base::DataVector(dimensions);
  EigenHelper::svd(EigenHelper::toEigen(matrix), i.eigenVectors, i.eigenValues);
  return i;
}

AsMcResult AsMcIntervalCutter::cut(const AsMcInput& input, const AsInfo& info) {
  size_t dimensions = info.eigenValues.size();
  std::mt19937_64 prng;
  std::uniform_int_distribution<size_t> dist(0, input.samples.getSize() - 1);

  std::vector<std::pair<double, double>> eigenValueIntervals(dimensions);
  for (size_t d = 0; d < dimensions; ++d) {
    eigenValueIntervals[d].first = info.eigenValues[d];
    eigenValueIntervals[d].second = info.eigenValues[d];
  }

  for (size_t i = 0; i < bootstrapSamples; i++) {
    sgpp::base::DataMatrix bootstrapMatrix(info.eigenValues.size(), info.eigenValues.size());
    for (size_t j = 0; j < input.samples.getSize(); i++) {
      size_t l = dist(prng);
      bootstrapMatrix.add(input.samples.getValues()[l]);
    }
    bootstrapMatrix.mult(1.0 / static_cast<double>(bootstrapSamples));

    sgpp::base::DataMatrix bootstrapEigenVectorMatrix(dimensions, dimensions);
    sgpp::base::DataVector bootstrapEigenValues(dimensions);
    EigenHelper::svd(EigenHelper::toEigen(bootstrapMatrix), bootstrapEigenVectorMatrix,
                     bootstrapEigenValues);
    for (size_t d = 0; d < dimensions; ++d) {
      double e = bootstrapEigenValues[d];
      if (e < eigenValueIntervals[d].first) {
        eigenValueIntervals[d].first = e;
      }
      if (e > eigenValueIntervals[d].second) {
        eigenValueIntervals[d].second = e;
      }
    }
  }

  double max = 0;
  size_t cutoff = 0;
  for (size_t d = 0; d < dimensions; ++d) {
    double size = eigenValueIntervals[d].second - eigenValueIntervals[d].first;
    if (size > max) {
      max = size;
      cutoff = d + 1;
    }
  }
  return AsMcResult(input, info.eigenVectors, cutoff, type, level, mean);
}

AsMcCutter::AsMcCutter(GridType type, level_t level, const DataVector& mean)
    : type(type), level(level), mean(mean) {}
}  // namespace base
}  // namespace sgpp
