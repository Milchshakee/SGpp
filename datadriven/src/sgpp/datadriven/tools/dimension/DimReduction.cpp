#include <eigen3/Eigen/Dense>
#include <sgpp/base/function/scalar/ChainScalarFunction.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/base/tools/EigenHelper.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/datadriven/tools/DatasetTools.hpp>
#include <sgpp/datadriven/tools/dimension/DimReduction.hpp>
#include <sgpp/base/function/scalar/SumFunction.hpp>

namespace sgpp {
namespace base {

double DimReduction::calculateMcL2Error(ScalarFunction& func, DistributionSample& dist) {
  static WrapperScalarFunction::FunctionEvalType zeroFunc = [](const DataVector&) { return 0; };
  WrapperScalarFunction zero(func.getNumberOfParameters(), zeroFunc);

  static WrapperVectorFunction::FunctionEvalType tFunc = [](const DataVector& v, DataVector& out) {
    out = v;
  };
  WrapperVectorFunction t(func.getNumberOfParameters(), func.getNumberOfParameters(), tFunc);
  return DimReduction::calculateMcL2Error(func, t, zero, dist);
}

PointSample<double> DimReduction::createActiveSubspaceSample(PointSample<double> input,
                                                             const DataMatrix& basis,
                                                             size_t reducedDims) {
  InputProjectionFunction f(basis, reducedDims);
  std::vector<DataVector> newPoints(input.getSize());
  for (size_t i = 0; i < input.getSize(); i++) {
    f.eval(input.getKeys()[i], newPoints[i]);
  }
  return PointSample<double>(newPoints, input.getValues());
}

datadriven::RegressionLearner getLearner(
    size_t dimension, sgpp::datadriven::RegularizationConfiguration regularizationConfig,
    DimReduction::RegressionConfig config) {
  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = dimension;
  gridConfig.level_ = config.gridLevel;

  gridConfig.type_ = sgpp::base::GridType::LinearBoundary;

  auto adaptivityConfig = sgpp::base::AdaptivityConfiguration();
  adaptivityConfig.numRefinementPoints_ = config.refinementPoints;
  adaptivityConfig.numRefinements_ = config.refinements;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = config.maxIterations;
  solverConfig.eps_ = 1e-11;
  solverConfig.threshold_ = 1e-8;
  solverConfig.verbose_ = true;

  return sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                             solverConfig, regularizationConfig);
}

std::string showRegularizationConfiguration(
    const sgpp::datadriven::RegularizationConfiguration& regularizationConfig) {
  std::ostringstream ss;
  const auto regType = regularizationConfig.type_;
  if (regType == sgpp::datadriven::RegularizationType::Diagonal) {
    ss << "type: DiagonalMatrix\t";
  } else if (regType == sgpp::datadriven::RegularizationType::Identity) {
    ss << "type: IdentityMatrix\t";
  } else if (regType == sgpp::datadriven::RegularizationType::Laplace) {
    ss << "type: Laplace\t";
  } else {
    ss << "type: unknown\t";
  }

  ss << "lambda: " << regularizationConfig.lambda_
     << "\tmultiplicationFactor: " << regularizationConfig.exponentBase_;
  return ss.str();
}

  void trainValidateSplit(datadriven::Dataset& data, double trainDataShare,
                        DataMatrix& trainX, DataVector& trainY,
                        DataMatrix& validateX,
                        DataVector& validateY)
{

      

std::vector<int> indices;
    for (int i = 0; i < data.getNumberInstances(); i++) indices.push_back(i);
  std::random_shuffle(indices.begin(), indices.end());

  size_t trainSize = trainDataShare * data.getTargets().size();
  trainX = DataMatrix(trainSize, data.getDimension());
  trainY = DataVector(trainSize);
  for (size_t i = 0; i < trainSize; i++) {
    DataVector row(data.getDimension());
    data.getData().getRow(indices[i], row);
    trainX.setRow(i, row);
    trainY[i] = data.getTargets()[indices[i]];
  }

  size_t valSize = data.getTargets().size() - trainSize;
  validateX = DataMatrix(valSize, data.getDimension());
  validateY = DataVector(valSize);
  for (size_t i = 0; i < valSize; i++) {
    DataVector row(data.getDimension());
    data.getData().getRow(indices[trainSize + i], row);
    validateX.setRow(i, row);
    validateY[i] = data.getTargets()[indices[trainSize + i]];
  }
}

void createRegressionSurrogate(datadriven::Dataset& data,
                                      DimReduction::RegressionConfig config, double totalError, SGridSample& out, size_t& gridPoints, double& error) {
  size_t trainSize = config.trainDataShare * data.getTargets().size();
  DataMatrix trainX = data.getData();
  trainX.resizeRows(trainSize);
  DataVector trainY = data.getTargets();
  trainY.resize(trainSize);

  size_t valSize = data.getTargets().size() - trainSize;
  DataMatrix validateX(valSize, data.getDimension());
  DataVector validateY(valSize);
  for (size_t i = 0; i < valSize; i++) {
    validateX[i] = data.getData()[trainSize + i];
    validateY[i] = data.getTargets()[trainSize + i];
  }

  datadriven::RegularizationConfiguration bestConfig;
  double bestError = std::numeric_limits<double>::infinity();
  for (double& lambda : config.lambdas) {
    for (double& expBase : config.regularizationBases) {
      const auto regularizationType = sgpp::datadriven::RegularizationType::Diagonal;
      auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
      regularizationConfig.type_ = regularizationType;
      regularizationConfig.lambda_ = lambda;
      regularizationConfig.exponentBase_ = expBase;
      regularizationConfig.regularizationMetric_ = datadriven::RegularizationMetricType::mse;

      double meanRMSE = 0;
      for (size_t i = 0; i < 10; i++) {
        auto learner = getLearner(data.getDimension(), regularizationConfig, config);
        learner.train(trainX, trainY);
        trainValidateSplit(data, config.trainDataShare, trainX, trainY, validateX, validateY);
        double curRMSE = std::sqrt(learner.getMSE(validateX, validateY));
        meanRMSE += curRMSE;
      }
      meanRMSE /= 10;
      double relError = meanRMSE / totalError;
      std::cout << "Tested parameters are\n" << showRegularizationConfiguration(regularizationConfig) << ".\n";
      if (relError < bestError) {
        std::cout << "Better! Relative L2 error is now " << relError << std::endl;
        bestError = relError;
        bestConfig = regularizationConfig;
      } else {
        std::cout << "Worse! Relative L2 error is now " << relError << std::endl;
      }
    }
  }

  error = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < 10; i++) {
    trainValidateSplit(data, config.trainDataShare * 4, trainX, trainY, validateX, validateY);
    auto bestLearner = getLearner(data.getDimension(), bestConfig, config);
      bestLearner.train(trainX, trainY);
      auto ptr = bestLearner.getGridPtr();
      SGridSample sample(ptr, bestLearner.getWeights());
    sample.setHierarchised(true);
    double curRMSE = std::sqrt(bestLearner.getMSE(validateX, validateY));
    double relError = curRMSE / totalError;
    if (relError < error) {
      std::cout << "Better! Relative L2 error is now " << relError << std::endl;
      gridPoints = sample.getSize();
      error = relError;
      out = sample;
    } 
  }
}

size_t asIntervalDimensions(DataVector& eigenValues, size_t bootstrapSamples,
                            size_t bootstrapCount,
                            PointSample<DataMatrix>& samples) {
  size_t dimensions = eigenValues.size();
  std::mt19937_64 prng;
  std::uniform_int_distribution<size_t> dist(0, samples.getSize() - 1);

  std::vector<std::pair<double, double>> eigenValueIntervals(dimensions);
  for (size_t d = 0; d < dimensions; ++d) {
    eigenValueIntervals[d].first = eigenValues[d];
    eigenValueIntervals[d].second = eigenValues[d];
  }

  for (size_t i = 0; i < bootstrapCount; i++) {
    sgpp::base::DataMatrix bootstrapMatrix(eigenValues.size(), eigenValues.size());
    for (size_t j = 0; j < bootstrapSamples; j++) {
      size_t l = dist(prng);
      bootstrapMatrix.add(samples.getValues()[l]);
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

  DataVector meanEigenValues(eigenValues.size());
  for (size_t d = 0; d < dimensions; ++d) {
    meanEigenValues[d] = (eigenValueIntervals[d].first + eigenValueIntervals[d].second) / 2.0;
  }

  double sum = meanEigenValues.sum();
  double partialSum = 0;
  double max = 0;
  size_t cutoff = 0;
  for (size_t d = 0; d < dimensions - 1; ++d) {
    partialSum += meanEigenValues[d];
    if (partialSum / sum > 0.9) {
      return d + 1;
      }
    //double size = (meanEigenValues[d] - meanEigenValues[d+1]) / meanEigenValues[d];
    //if (size > max) {
    //  max = size;
    //  cutoff = d + 1;
    //}
  }
  return cutoff;
}

PointSample<DataMatrix> fromFiniteDifferences(ScalarFunction& func, DistributionSample& v,
                                              double h) {
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

ActiveSubspaceInfo activeSubspaceMC(
    sgpp::base::PointSample<sgpp::base::DataMatrix>& m) {
  size_t dimensions = m.getDimensions();
  sgpp::base::DataMatrix matrix(dimensions, dimensions);
  for (size_t i = 0; i < m.getSize(); ++i) {
    matrix.add(m.getValues()[i]);
  }
  matrix.mult(1.0 / static_cast<double>(m.getSize()));

  ActiveSubspaceInfo i;
  i.eigenVectors = sgpp::base::DataMatrix(dimensions, dimensions);
  i.eigenValues = sgpp::base::DataVector(dimensions);
  EigenHelper::svd(EigenHelper::toEigen(matrix), i.eigenVectors, i.eigenValues);
  return i;
}


AsReductionResult DimReduction::reduceAS(std::shared_ptr<ScalarFunction>& f,
                                       sgpp::base::DistributionsVector dist,
                                       std::vector<RegressionConfig> configs) {
  std::vector<GridReductionResult> results;
  std::shared_ptr<ScalarFunction> current = f;
  std::vector<std::shared_ptr<ScalarFunction>> funcs;
  size_t counter = 0;
  double currentError = 1;
  size_t gridPoints = 0;
  auto errorSample = sgpp::base::DistributionSample(configs[0].errorCalcSamples, dist);
  double totalError = calculateMcL2Error(*f, errorSample);
  for (RegressionConfig& config : configs) {
    std::cout << std::endl;
    std::cout << "Starting iteration: " << counter << std::endl;

    auto distSample = sgpp::base::DistributionSample(config.samples, dist);
    sgpp::base::PointSample<sgpp::base::DataMatrix> m = fromFiniteDifferences(*current, distSample, 0.00000001);
    auto i = activeSubspaceMC(m);
    size_t reducedDims = 2;
    //asIntervalDimensions(i.eigenValues, 0.5 * distSample.getSize(), 5, m);
    i.eigenVectors.resizeRowsCols(i.eigenVectors.getNrows(), reducedDims);

    std::cout << "Eigenvalues " << i.eigenValues.toString() << "\n";
    std::cout << "Reducing to " << reducedDims << " dimensions.\n\n";

    sgpp::base::PointSample<double> sample = sgpp::base::DimReduction::createActiveSubspaceSample(
        sgpp::base::SampleHelper::sampleScalarFunction(distSample, *current),
        i.eigenVectors, reducedDims);
    sgpp::datadriven::Dataset data = sgpp::base::SampleHelper::fromPointSample(sample);

    sgpp::base::GridReductionResult result = sgpp::base::DimReduction::activeSubspaceReduction(
        current, data, i.eigenVectors, reducedDims, config, totalError);
    results.push_back(result);
    funcs.push_back(result.replacementFunction);

    std::shared_ptr<ScalarFunction> sum = std::make_shared<SumFunction>(funcs);
    auto fullResult = AsReductionResult{f, sum, results};

    gridPoints += result.reducedSample.getSize();
    currentError = (fullResult.mcL2Error(errorSample) / totalError);
    std::cout << "New relative L2 error: \n"
              << currentError << ".\n";
    std::cout << "Grid points: " << gridPoints << ".\n";
    current = result.errorFunction;
    counter++;
  }
  std::shared_ptr<ScalarFunction> sum = std::make_shared<SumFunction>(funcs);
  return AsReductionResult{f, sum, results};
}

ReductionResult::ReductionResult(std::shared_ptr<ScalarFunction>& func,
                                 std::shared_ptr<ScalarFunction> replacement)
    : originalFunction(func), replacementFunction(replacement) {
  std::vector<std::shared_ptr<ScalarFunction>> fs{func, replacement};
  errorFunction = std::make_shared<SumFunction>(fs, std::vector<bool>{true, false});
}


double ReductionResult::mcL2Error(DistributionSample& sample) {
  double res = 0;
  for (size_t i = 0; i < sample.getSize(); i++) {
    double val = originalFunction->eval(sample.getVectors()[i]);
    double newVal = replacementFunction->eval(sample.getVectors()[i]);
    res += pow(val - newVal, 2);
  }

  return sqrt(res / static_cast<double>(sample.getSize()));
}

GridReductionResult::GridReductionResult(std::shared_ptr<ScalarFunction>& func,
                                         std::shared_ptr<VectorFunction>& t,
                                         std::shared_ptr<ScalarFunction>& sampleFunction,
                                         SGridSample& sample, size_t gridPoints, double l2Error)
    : ReductionResult(func, std::make_shared<ChainScalarFunction>(
                                std::vector<std::shared_ptr<VectorFunction>>(1, t),
                                sampleFunction)),
      reducedSample(sample),
      reducedFunction(sampleFunction),
      transformation(t), gridPoints(gridPoints), l2Error(l2Error) {
}


AsReductionResult::AsReductionResult(
    std::shared_ptr<ScalarFunction>& func,
    std::shared_ptr<ScalarFunction> replacement, std::vector<GridReductionResult>& results) : ReductionResult(func, replacement), reductions(results){}

GridReductionResult DimReduction::activeSubspaceReduction(std::shared_ptr<ScalarFunction>& f,
                                                      datadriven::Dataset& data,
                                                      const DataMatrix& basis, size_t reducedDims,
                                                      RegressionConfig config, double totalError) {
  SGridSample s;
  size_t gridPoints;
  double error;
  createRegressionSurrogate(data, config, totalError, s, gridPoints, error);
  std::shared_ptr<VectorFunction> projection =
      std::make_shared<InputProjectionFunction>(basis, reducedDims);
  std::shared_ptr<ScalarFunction> reducedFunc = std::make_shared<InterpolantScalarFunction>(s);
  return GridReductionResult{f, projection, reducedFunc, s, gridPoints, error};
}

double sgpp::base::DimReduction::calculateMcL2Error(ScalarFunction& func,
                                                    VectorFunction& transformation,
                                                    ScalarFunction& reduced,
                                                    DistributionSample& dist) {
  size_t funcDimensions = func.getNumberOfParameters();
  size_t newDimensions = reduced.getNumberOfParameters();

  sgpp::base::DataVector point(funcDimensions);
  double res = 0;

  for (size_t i = 0; i < dist.getSize(); i++) {
    point = dist.getVectors()[i];
    double val = func.eval(point);
    DataVector out(newDimensions);
    transformation.eval(point, out);
    double newVal = reduced.eval(out);
    res += pow(val - newVal, 2);
  }

  return sqrt(res / static_cast<double>(dist.getSize()));
}


std::shared_ptr<InputProjectionFunction> InputProjectionFunction::identity(size_t dims) {
  sgpp::base::DataMatrix id(dims, dims);
  for (size_t d = 0; d < dims; ++d) {
    sgpp::base::DataVector unit(dims, 0.0);
    unit[d] = 1;
    id.setColumn(d, unit);
  }
  return std::make_shared<InputProjectionFunction>(id, dims);
}

InputProjectionFunction::InputProjectionFunction(const DataMatrix& as, size_t reducedDims)
    : VectorFunction(as.getNrows(), reducedDims),
      activeSubspaceMat(as),
      localBasis(reducedDims, reducedDims) {
  size_t dims = as.getNrows();
  activeSubspaceMat.transpose();
  DataMatrix zonotope(reducedDims, dims);

  for (size_t d = 0; d < dims; ++d) {
    DataVector unit(dims, 0.0);
    unit[d] = 1;
    DataVector generator = EigenHelper::mult(activeSubspaceMat, unit);
    zonotope.setColumn(d, generator);
  }

  auto eigenBasis = EigenHelper::toEigen(zonotope);
  Eigen::MatrixXd q = eigenBasis.householderQr().householderQ();
  DataMatrix ortho = EigenHelper::fromEigen(q);

  for (size_t d = 0; d < reducedDims; ++d) {
    DataVector col(reducedDims);
    ortho.getColumn(d, col);

    double sum = 0;
    for (size_t k = 0; k < dims; ++k) {
      DataVector z(reducedDims);
      zonotope.getColumn(k, z);
      sum += std::abs(col.dotProduct(z));
    }
    col.mult(1 / sum);
    localBasis.setColumn(d, col);
  }
  localBasis.transpose();
}

InputProjectionFunction::~InputProjectionFunction() {}

void InputProjectionFunction::eval(const DataVector& x, DataVector& value) {
  DataVector v = x;
  v.sub(DataVector(x.getSize(), 0.5));
  value = EigenHelper::mult(activeSubspaceMat, v);
  value = EigenHelper::mult(localBasis, value);
  value.add(DataVector(value.getSize(), 0.5));
  for (size_t d = 0; d < value.size(); ++d) {
    if (value[d] < -1e-8 || value[d] > 1 + 1e-8) {
      throw std::logic_error("Out of bounds");
    }

    if (value[d] < 0) {
      value[d] = 0;
    }
    if (value[d] > 1) {
      value[d] = 1;
    }
  }
  }

void InputProjectionFunction::clone(std::unique_ptr<VectorFunction>& clone) const {}

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

ReductionResult DimReduction::reduce(ScalarFunction& f, const DataMatrix& basis, size_t reducedDims,
                                     AnovaTypes::level_t level) {
  WrapperScalarFunction::FunctionEvalType func = [&reducedDims, &basis, &f](const DataVector& v) {
    return 0;
  };
  WrapperScalarFunction wrapped(f.getNumberOfParameters(), func);

  std::shared_ptr<sgpp::base::Grid> grid(
      sgpp::base::Grid::createAnovaPrewaveletBoundaryGrid(f.getNumberOfParameters()));
  grid->getGenerator().regular(AnovaTypes::toNormalLevel(level));
  sgpp::base::SGridSample sample(grid, wrapped);
  sample.hierarchise();
  SGridSample reducedSample = createReducedAnovaSample(sample, level, reducedDims);
  auto projection = std::make_shared<InputProjectionFunction>(basis, reducedDims);
  // return {f, projection, reducedSample};
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

}  // namespace base
}  // namespace sgpp
