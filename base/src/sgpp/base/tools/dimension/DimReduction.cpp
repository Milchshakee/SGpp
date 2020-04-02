#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/tools/EigenHelper.hpp>
#include <sgpp/base/tools/dimension/AsMcReducer.hpp>
#include <sgpp/base/tools/dimension/DimReduction.hpp>
#include <sgpp/base/tools/dimension/HypercubeFitFunction.hpp>
#include <sgpp/base/tools/dist/RandomUniformDistribution.hpp>

namespace sgpp {
namespace base {

double sgpp::base::DimReduction::calculateMcL2Error(ScalarFunction& func,
                                                    VectorFunction& transformation,
                                                    ScalarFunction& reduced, uint64_t seed,
                                                    size_t samples) {
  std::mt19937_64 rand(seed);
  std::uniform_real_distribution<double> dist(0, 1);
  size_t funcDimensions = func.getNumberOfParameters();
  size_t newDimensions = reduced.getNumberOfParameters();

  sgpp::base::DataVector point(funcDimensions);
  double res = 0;

  for (size_t i = 0; i < samples; i++) {
    for (size_t d = 0; d < funcDimensions; d++) {
      point[d] = dist(rand);
    }
    double val = func.eval(point);
    DataVector out(newDimensions);
    transformation.eval(point, out);
    double newVal = reduced.eval(out);
    res += pow(val - newVal, 2);
  }

  return sqrt(res / static_cast<double>(samples));
}

InputProjectionFunction::InputProjectionFunction(const DataMatrix& basis, size_t reducedDims)
    : VectorFunction(basis.getNcols(), reducedDims) {
  DataMatrix copy = basis;
  copy.transpose();
  copy.resize(reducedDims);
  InputProjectionFunction::basis = copy;

  size_t dims = basis.getNcols();
    posRanges = DataVector(reducedDims, 0);
  double m = 0.5;
  for (size_t d = 0; d < reducedDims; ++d) {
    DataVector col(dims);
    basis.getColumn(d, col);
    for (size_t j = 0; j < dims; ++j) {
      double v = col[j];

      posRanges[d] += v >= 0 ? (1 - m) * v : -m * v;
    }
  }
}

InputProjectionFunction::~InputProjectionFunction() {}

void InputProjectionFunction::eval(const DataVector& x, DataVector& value) {
  DataVector v = x;
  v.sub(DataVector(x.getSize(), 0.5));
  value = EigenHelper::mult(basis, v);
  for (size_t d = 0; d < basis.getNrows(); ++d) {
    value[d] *= posRanges[d];
  }
  value.add(DataVector(basis.getNrows(), 0.5));
}

void InputProjectionFunction::clone(std::unique_ptr<VectorFunction>& clone) const {}

DataVector DimReduction::transformPoint(const DataMatrix& basis, const DataVector& in,
                                        size_t reducedDims) {
  size_t dims = in.getSize();
  DataVector point = in;

  DataVector posRanges(reducedDims, 0);
  double m = 0.5;
  for (size_t d = 0; d < reducedDims; ++d) {
    DataVector col(dims);
    basis.getColumn(d, col);
    for (size_t j = 0; j < dims; ++j) {
      double v = col[j];

      posRanges[d] += v >= 0 ? (1 - m) * v : -m * v;
    }
  }
  for (size_t d = 0; d < reducedDims; ++d) {
    point[d] = m + ((in[d] - m) * (2 * posRanges[d]));
  }

  DataVector relative = point;
  DataVector base = point;
  for (size_t d = reducedDims; d < dims; ++d) {
    base.set(d, 0.5);
  }
  relative.sub(base);
  relative.mult(1.0 / std::sqrt(0.5 * 0.5 * (dims - reducedDims)));

  base.sub(DataVector(dims, 0.5));
  base = EigenHelper::mult(basis, base);
  base.add(DataVector(dims, 0.5));
  relative = EigenHelper::mult(basis, relative);

  if (relative.getNumberNonZero() > 0) {
    double toScale = calcScalingFactor(base, relative);
    relative.mult(toScale);
    base.add(relative);
  }

  // cap
  for (size_t d = 0; d < dims; ++d) {
    if (base[d] < 0) {
      base[d] = 0;
    } else if (base[d] > 1) {
      base[d] = 1;
    }
  }

  return base;
}

double DimReduction::calcScalingFactor(sgpp::base::DataVector& point,
                                       sgpp::base::DataVector& direction) {
  size_t dim = point.getSize();
  double current = std::numeric_limits<double>::infinity();
  for (size_t d = 0; d < dim; ++d) {
    double zero = -point[d] / direction[d];
    double one = (1 - point[d]) / direction[d];
    double max = std::max(zero, one);
    if (max != -std::numeric_limits<double>::infinity() && max < current) {
      current = max;
    }
  }
  return current;
}

sgpp::base::SGridSample DimReduction::createReducedAnovaSample(sgpp::base::SGridSample& sample,
                                                               AnovaTypes::level_t level,
                                                               size_t reducedDims) {
  std::shared_ptr<sgpp::base::Grid> grid(
      sgpp::base::Grid::createAnovaPrewaveletBoundaryGrid(reducedDims));
  grid->getGenerator().regular(AnovaTypes::toNormalLevel(level));
  std::function<double(const DataVector&)> func = [&sample](const DataVector& v) {
    DataVector coords = v;
    coords.resizeZero(sample.getDimensions());
    double value = sample.getValue(coords);
    return value;
  };

  SGridSample reducedSample(grid, func);
  reducedSample.setHierarchised(true);
  return reducedSample;
}

ActiveSubspaceInfo DimReduction::activeSubspaceMC(ScalarFunction& f, VectorDistribution& dist) {
  sgpp::base::PointSample<sgpp::base::DataMatrix> m =
      sgpp::base::AsMcReducer::fromFiniteDifferences(f, dist, 0.0001);
  size_t dimensions = f.getNumberOfParameters();
  sgpp::base::DataMatrix matrix(dimensions, dimensions);
  for (size_t i = 0; i < dist.getSize(); ++i) {
    matrix.add(m.getValues()[i]);
  }
  matrix.mult(1.0 / static_cast<double>(dist.getSize()));

  ActiveSubspaceInfo i;
  i.eigenVectors = sgpp::base::DataMatrix(dimensions, dimensions);
  i.eigenValues = sgpp::base::DataVector(dimensions);
  EigenHelper::svd(EigenHelper::toEigen(matrix), i.eigenVectors, i.eigenValues);
  return i;
}

ReductionResult DimReduction::reduce(ScalarFunction& f, const DataMatrix& basis, size_t reducedDims,
                                     AnovaTypes::level_t level) {
  WrapperScalarFunction::FunctionEvalType func = [&reducedDims, &basis, &f](const DataVector& v) {
    DataVector transformed = transformPoint(basis, v, reducedDims);
    double value = f.eval(transformed);
    if (std::isnan(value)) {
      value = 0;
    }
    return value;
  };
  WrapperScalarFunction wrapped(f.getNumberOfParameters(), func);

  std::shared_ptr<sgpp::base::Grid> grid(
      sgpp::base::Grid::createAnovaPrewaveletBoundaryGrid(f.getNumberOfParameters()));
  grid->getGenerator().regular(AnovaTypes::toNormalLevel(level));
  sgpp::base::SGridSample sample(grid, wrapped);
  sample.hierarchise();
  SGridSample reducedSample = createReducedAnovaSample(sample, level, reducedDims);
  return {InputProjectionFunction(basis, reducedDims), reducedSample};
}

InputProjection::InputProjection(const DataMatrix& basis, size_t n, const DataVector& mean)
    : oldDimensions(basis.getNrows()), newDimensions(n), mean(mean), func(*this) {
  oldToNewBasis = basis;
  oldToNewBasis.transpose();
  oldToNewBasis.resizeRowsCols(newDimensions, oldDimensions);
  newToOldBasis = basis;
  newToOldBasis.transpose();
  newToOldBasis.resizeRowsCols(oldDimensions, newDimensions);
  calculateRanges();
}

void InputProjection::calculateRanges() {
  size_t dim = newDimensions;
  posRange = DataVector(dim, 0);
  negRange = DataVector(dim, 0);
  for (size_t d = 0; d < dim; ++d) {
    for (size_t j = 0; j < oldDimensions; ++j) {
      DataVector col(oldDimensions);
      oldToNewBasis.getRow(d, col);
      double v = col[j];
      double m = mean[j];

      posRange[d] += v >= 0 ? (1 - m) * v : -m * v;
      negRange[d] += v < 0 ? (1 - m) * v : -m * v;
    }
  }

  start = mean;
  end = mean;
  for (size_t d = 0; d < newDimensions; ++d) {
    DataVector add(oldDimensions);
    oldToNewBasis.getRow(d, add);
    add.mult(negRange[d]);
    start.add(add);
    oldToNewBasis.getRow(d, add);
    add.mult(posRange[d]);
    end.add(add);
  }
}

void InputProjection::inverse(const DataVector& in, DataVector& out) {
  out = start;
  for (size_t d = 0; d < newDimensions; ++d) {
    double scale = posRange[d] - negRange[d];
    DataVector add(oldDimensions);
    oldToNewBasis.getRow(d, add);
    add.mult(scale * in[d]);
    out.add(add);
  }
}

size_t InputProjection::getNewDimensions() { return newDimensions; }

const DataMatrix InputProjection::getTransformationMatrix() { return oldToNewBasis; }

const DataVector& InputProjection::getStart() { return start; }

const DataVector& InputProjection::getEnd() { return end; }

InputProjection::ProjectionFunction::ProjectionFunction(InputProjection& p)
    : VectorFunction(p.oldDimensions, p.newDimensions), p(&p) {}

void InputProjection::ProjectionFunction::eval(const DataVector& in, DataVector& out) {
  out = in;
  out.sub(p->mean);
  out = EigenHelper::mult(p->oldToNewBasis, out);
  DataVector ranges(p->newDimensions);
  for (size_t d = 0; d < p->newDimensions; ++d) {
    double scale = p->posRange[d] - p->negRange[d];
    if (out[d] > p->posRange[d]) {
      out[d] = p->posRange[d];
    } else if (out[d] < p->negRange[d]) {
      out[d] = p->negRange[d];
    }

    out[d] = (out[d] - p->negRange[d]) / scale;
  }
}

void InputProjection::ProjectionFunction::clone(std::unique_ptr<VectorFunction>& clone) const {}

VectorFunction& InputProjection::getFunction() { return func; }

L2SquaredMcRule::L2SquaredMcRule(uint64_t seed, size_t samples) : seed(seed), samples(samples) {}

L2McRule::L2McRule(uint64_t seed, size_t samples) : seed(seed), samples(samples) {}

double L2McRule::calculateRelativeError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r) {
  double e = calculateAbsoluteError(f, t, r);
  double fe = calculateAbsoluteError(f);
  return fe != 0.0 ? e / fe : 0.0;
}

double L2McRule::calculateAbsoluteError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r) {
  OperationL2 o(seed, samples);
  double e = o.calculateMcL2Error(f, t, r);
  return e;
}

double L2McRule::calculateAbsoluteError(ScalarFunction& f) {
  OperationL2 o(seed, samples);
  double l2 = o.calculateMcL2Norm(f);
  return l2;
}

double ErrorRule::calculateRelativeError(ScalarFunction& f, VectorFunction& t, ScalarFunction& r) {
  double e = calculateAbsoluteError(f, t, r);
  double fe = calculateAbsoluteError(f);
  return fe != 0.0 ? e / fe : 0.0;
}

double L2SquaredMcRule::calculateAbsoluteError(ScalarFunction& f, VectorFunction& t,
                                               ScalarFunction& r) {
  OperationL2 o(seed, samples);
  double e = o.calculateMcL2Error(f, t, r);
  return std::pow(e, 2);
}

double L2SquaredMcRule::calculateAbsoluteError(ScalarFunction& f) {
  OperationL2 o(seed, samples);
  double l2 = o.calculateMcL2Norm(f);
  return std::pow(l2, 2);
}
}  // namespace base
}  // namespace sgpp
