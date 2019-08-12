

#include <sgpp/base/tools/dimension/AsMcReducer.hpp>

namespace sgpp {
namespace base {


ScalarFunction& AsMcResult::getOriginalFunction() { return *originalFunc; }

VectorFunction& AsMcResult::getTransformationFunction() { return projection.getFunction(); }

AsMcFixedCutter::AsMcFixedCutter(size_t n, GridType type, level_t level)
    : FixedCutter<sgpp::base::AsMcInput, sgpp::base::AsInfo, sgpp::base::AsMcResult>(n), AsMcCutter(type, level) {}

AsMcErrorRuleCutter::AsMcErrorRuleCutter(ErrorRule& r, double maxError, GridType type,
                                         level_t level)
    : ErrorRuleCutter<sgpp::base::AsMcInput, sgpp::base::AsInfo, sgpp::base::AsMcResult>(
          r, maxError), AsMcCutter(type, level) {}

AsMcResult AsMcErrorRuleCutter::cut(const AsMcInput& input, const AsInfo& info) {
  AsMcResult last(input, info.eigenVectors, input.function.getNumberOfParameters(), type, level);
  for (size_t d = 1; d < input.function.getNumberOfParameters(); ++d) {
    AsMcResult result(input, info.eigenVectors, input.function.getNumberOfParameters() - d, type, level);
    if (result.calculateRelativeError(r) > maxError) {
      return last;
    } else {
      last = result;
    }
  }
  return last;
}

ScalarFunction& AsMcResult::getReducedFunction() { return evalFunc; }

SGridSample& AsMcResult::getReducedOutput() { return reduced; }

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
  VectorDistribution& v) {
}


AsMcResult::AsMcResult(const AsMcInput& input, const DataMatrix& m, size_t n, GridType type,
                       level_t l)
    : projection(m, n, DataVector(m.getNrows(), 0.5)) {
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
  originalFunc = &input.function;
}

AsMcResult AsMcFixedCutter::cut(const AsMcInput& input, const AsInfo& info) {
  return AsMcResult(input, info.eigenVectors, n, type, level);
}

AsMcIntervalCutter::AsMcIntervalCutter(size_t bootstrapSamples, GridType type, level_t level)
    : AsMcCutter(type, level), bootstrapSamples(bootstrapSamples) {}

  
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
    EigenHelper::svd(EigenHelper::toEigen(bootstrapMatrix), bootstrapEigenVectorMatrix, bootstrapEigenValues);
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
  return AsMcResult(input, info.eigenVectors, cutoff, type, level);
}


AsMcCutter::AsMcCutter(GridType type, level_t level) : type(type), level(level) {
}
}  // namespace base
}  // namespace sgpp
