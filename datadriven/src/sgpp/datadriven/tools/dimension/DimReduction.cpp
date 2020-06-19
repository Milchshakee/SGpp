#include <eigen3/Eigen/Dense>
#include <sgpp/base/function/scalar/ChainScalarFunction.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/SumScalarFunction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/function/vector/SumVectorFunction.hpp>
#include <sgpp/base/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/base/tools/EigenHelper.hpp>
#include <sgpp/base/tools/sle/solver/Eigen.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrixGradientMC.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/datadriven/tools/DatasetTools.hpp>
#include <sgpp/datadriven/tools/dimension/DimReduction.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>

namespace sgpp {
namespace base {

    
DimReduction::L2 DimReduction::L2_ERROR = {};
DimReduction::L2 DimReduction::RMSE_ERROR = {};
DimReduction::MSE DimReduction::MSE_ERROR = {};
DimReduction::NRMSE DimReduction::NRMSE_ERROR = {};


std::tuple<PointSample<double>, PointSample<double>> createActiveSubspaceSample(
    PointSample<double>& originalInput, PointSample<double>& input, const DataMatrix& basis,
    size_t sampleCount) {
  std::vector<int> indices(input.getSize());
  for (int j = 0; j < input.getSize(); j++) indices[j] = j;
  std::random_shuffle(indices.begin(), indices.end());

  InputProjectionFunction f(basis, basis.getNcols());
  std::vector<DataVector> newPoints(sampleCount);
  std::vector<double> newValues(sampleCount);
  for (size_t i = 0; i < sampleCount; i++) {
    f.eval(input.getKeys()[indices[i]], newPoints[i]);
    newValues[i] = input.getValues()[indices[i]];
  }

  std::vector<DataVector> testPoints(input.getSize() - sampleCount);
  std::vector<double> testValues(input.getSize() - sampleCount);
  for (size_t i = sampleCount; i < input.getSize(); i++) {
    testPoints[i - sampleCount] = originalInput.getKeys()[indices[i]];
    testValues[i - sampleCount] = originalInput.getValues()[indices[i]];
  }
  return std::tuple<PointSample<double>, PointSample<double>>{
      PointSample<double>(newPoints, newValues), PointSample<double>(testPoints, testValues)};
}

datadriven::RegressionLearner getLearner(
    size_t dimension, sgpp::datadriven::RegularizationConfiguration regularizationConfig,
    DimReduction::RegressionConfig config) {
  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = dimension;
  gridConfig.level_ = std::max<double>((std::log2(config.maxGridPoints) / dimension) - 1, 0);

  gridConfig.type_ = sgpp::base::GridType::LinearBoundary;

  auto adaptivityConfig = sgpp::base::AdaptivityConfiguration();
  adaptivityConfig.numRefinementPoints_ =
      config.refinementPointsPerGridPoint * config.maxGridPoints;
  adaptivityConfig.numRefinements_ = config.refinements;
  adaptivityConfig.numCoarseningPoints_ = 0;
  adaptivityConfig.coarsenInitialPoints_ = false;
  // adaptivityConfig.errorBasedRefinement = true;
  // adaptivityConfig.percent_ = 1000;
  // adaptivityConfig.refinementPeriod = 5;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::BiCGSTAB;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-11;
  solverConfig.threshold_ = 1e-8;
  solverConfig.verbose_ = true;

  return sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                             solverConfig, regularizationConfig);
}

void trainValidateSplit(PointSample<double>& data, size_t trainData, DataMatrix& trainX,
                        DataVector& trainY, DataMatrix& validateX, DataVector& validateY) {
  std::vector<int> indices;
  for (int i = 0; i < data.getSize(); i++) indices.push_back(i);
  std::random_shuffle(indices.begin(), indices.end());

  size_t trainSize = trainData;
  trainX = DataMatrix(trainSize, data.getDimensions());
  trainY = DataVector(trainSize);
  for (size_t i = 0; i < trainSize; i++) {
    trainX.setRow(i, data.getKeys()[indices[i]]);
    trainY[i] = data.getValues()[indices[i]];
  }

  size_t valSize = data.getSize() - trainSize;
  validateX = DataMatrix(valSize, data.getDimensions());
  validateY = DataVector(valSize);
  for (size_t i = 0; i < valSize; i++) {
    validateX.setRow(i, data.getKeys()[indices[trainSize + i]]);
    validateY[i] = data.getValues()[indices[trainSize + i]];
  }
}

double crossValidate(PointSample<double>& data, DataMatrix& trainX, DataVector& trainY,
                     DataMatrix& validateX, DataVector& validateY,
                     DimReduction::RegressionConfig& config,
                     datadriven::RegularizationConfiguration& regularizationConfig) {
  double meanRMSE = 0;
  for (size_t i = 0; i < config.crossValidations; i++) {
    auto learner = getLearner(data.getDimensions(), regularizationConfig, config);
    learner.train(trainX, trainY);
    trainValidateSplit(data, config.trainDataPerGridPoint * config.maxGridPoints, trainX, trainY,
                       validateX, validateY);
    double curRMSE = std::sqrt(learner.getMSE(validateX, validateY));
    meanRMSE += curRMSE;
  }
  meanRMSE /= config.crossValidations;
  return meanRMSE;
}

void createRegressionSurrogate(PointSample<double>& data, DimReduction::RegressionConfig config,
                               SGridSample& out, double& error) {
  DataMatrix trainX;
  DataVector trainY;
  DataMatrix validateX;
  DataVector validateY;

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

      double meanRMSE =
          crossValidate(data, trainX, trainY, validateX, validateY, config, regularizationConfig);
      // std::cout << "Tested parameters are\n" <<
      // showRegularizationConfiguration(regularizationConfig) << ".\n";
      if (meanRMSE < bestError) {
        // std::cout << "Better! Relative L2 error is now " << relError << std::endl;
        bestError = meanRMSE;
        bestConfig = regularizationConfig;
      } else {
        // std::cout << "Worse! Relative L2 error is now " << relError << std::endl;
      }
    }
  }

  error = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < 5; i++) {
    trainValidateSplit(data, config.trainDataPerGridPoint * config.maxGridPoints, trainX, trainY,
                       validateX, validateY);
    auto bestLearner = getLearner(data.getDimensions(), bestConfig, config);
    bestLearner.train(trainX, trainY);
    auto ptr = bestLearner.getGridPtr();
    double curRMSE = std::sqrt(bestLearner.getMSE(validateX, validateY));
    if (curRMSE < error) {
      // std::cout << "Better! Relative L2 error is now " << relError << std::endl;

      SGridSample sample(ptr, bestLearner.getWeights());
      sample.setHierarchised(true);
      error = curRMSE;
      out = sample;
    }
  }
}

size_t asIntervalDimensions(DataVector& eigenValues, size_t bootstrapSamples, size_t bootstrapCount,
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
    // double size = (meanEigenValues[d] - meanEigenValues[d+1]) / meanEigenValues[d];
    // if (size > max) {
    //  max = size;
    //  cutoff = d + 1;
    //}
  }
  return cutoff;
}

std::shared_ptr<VectorFunction> DimReduction::finiteDifferencesFunction(
    std::shared_ptr<ScalarFunction>& f, double h) {
  WrapperVectorFunction::FunctionEvalType func = [&f, h](const DataVector& vec, DataVector& out) {
    DataVector working = vec;
    for (size_t d = 0; d < f->getNumberOfParameters(); ++d) {
      double hOffset = vec[d] + h;
      if (hOffset > 1.0) {
        hOffset -= 2 * h;
      }
      working[d] = hOffset;
      double val = (f->eval(vec) - f->eval(working)) / h;
      out[d] = val;
      working[d] = vec[d];
    }
  };
  return std::make_shared<WrapperVectorFunction>(f->getNumberOfParameters(),
                                                 f->getNumberOfParameters(), func);
}

PointSample<DataVector> createGradientSample(PointSample<double>& v,
                                             DimReduction::GradientConfig& gradConfig) {
  std::vector<DataVector> keys(gradConfig.sampleCount);
  std::vector<DataVector> gradients(gradConfig.sampleCount);

  std::vector<int> indices(v.getSize());
  for (size_t j = 0; j < v.getSize(); j++) {
    indices[j] = j;
  }
  std::random_shuffle(indices.begin(), indices.end());
  for (size_t j = 0; j < gradConfig.sampleCount; j++) {
    keys[j] = v.getKeys()[indices[j]];
  }

  if (gradConfig.type == DimReduction::GradientConfig::GRADIENT_FUNCTION) {
    sgpp::base::DataVector sampleGradient(v.getDimensions());
    DataVector working(v.getDimensions());
    for (size_t i = 0; i < gradConfig.sampleCount; ++i) {
      const DataVector& vec = v.getKeys()[indices[i]];
      gradConfig.gradientFunction->eval(vec, working);
      gradients[i] = working;
    }
  } else if (gradConfig.type == DimReduction::GradientConfig::RANDOM_NEIGHBOUR_APPROXIMATION ||
             gradConfig.type == DimReduction::GradientConfig::NEAREST_NEIGHBOUR_APPROXIMATION) {
    for (size_t i = 0; i < gradConfig.sampleCount; ++i) {
      std::vector<DataVector> neighbours;
      std::vector<double> neighbourValues;

      for (size_t j = 0; j < gradConfig.neighbourConsiderationCount; ++j) {
        bool shouldSelect = true;
        if (gradConfig.type == DimReduction::GradientConfig::NEAREST_NEIGHBOUR_APPROXIMATION) {
          DataVector vec = v.getKeys()[i];
          vec.sub(v.getKeys()[indices[j]]);
          vec.abs();
          for (size_t d = 0; d < vec.getSize(); ++d) {
            vec[d] = std::pow(vec[d], gradConfig.pNorm);
          }
          double norm = std::pow(vec.sum(), 1.0 / gradConfig.pNorm);
          if (norm > gradConfig.maxNorm) {
            shouldSelect = false;
          }
        }

        if (shouldSelect) {
          neighbours.push_back(v.getKeys()[indices[j]]);
          neighbourValues.push_back(v.getValues()[indices[j]]);
        }
      }

      if (neighbours.size() == 0) continue;

      DataMatrix A(neighbours.size(), v.getDimensions());
      DataVector b(neighbours.size());
      for (size_t j = 0; j < neighbours.size(); ++j) {
        DataVector vec = v.getKeys()[i];
        vec.sub(neighbours[j]);
        A.setRow(j, vec);
        b[j] = neighbourValues[j] - v.getValues()[i];
      }

      EigenHelper::solveSLE(A, b, gradients[i]);
    }
  }
  return PointSample<DataVector>(keys, gradients);
}

DataMatrix randomBasis(size_t dims, size_t reducedDims) {
  DataMatrix mat(dims, reducedDims);
  mat.setAll(0.0);
  std::uniform_real_distribution<double> uniform(-1, 1);
  auto gen = std::default_random_engine();
  for (size_t i = 0; i < reducedDims; ++i) {
    DataVector randomVec(dims);
    for (size_t d = 0; d < dims; ++d) {
      randomVec[d] = uniform(gen);
    }
    randomVec.mult(1.0 / randomVec.l2Norm());
    mat.setColumn(i, randomVec);
  }

  auto eigenBasis = EigenHelper::toEigen(mat);
  Eigen::MatrixXd q = eigenBasis.householderQr().householderQ();
  DataMatrix ortho = EigenHelper::fromEigen(q);
  ortho.transpose();
  ortho.resizeRows(reducedDims);
  ortho.transpose();
  return ortho;
}

DataMatrix createBasis(PointSample<DataVector>& samples, DimReduction::BasisConfig& bConfig,
                       DimReduction::ReductionConfig& redConfig, DimReduction::Output& output) {
  size_t dimensions = samples.getDimensions();
  DataMatrix mat(samples.getDimensions(), samples.getDimensions());

  if (bConfig.type == DimReduction::BasisConfig::RANDOM) {
    return randomBasis(samples.getDimensions(), redConfig.maxReducedDimension);
  }

  if (bConfig.type == DimReduction::BasisConfig::ACTIVE_SUBSPACE) {
    for (size_t i = 0; i < samples.getSize(); ++i) {
      DataMatrix local(samples.getDimensions(), samples.getDimensions());
      for (size_t d = 0; d < samples.getDimensions(); ++d) {
        sgpp::base::DataVector col = samples.getValues()[i];
        col.mult(samples.getValues()[i][d]);
        local.setColumn(d, col);
      }
      mat.add(local);
    }
    mat.mult(1.0 / (samples.getSize()));
  } else if (bConfig.type == DimReduction::BasisConfig::INV_AS) {
    DataVector gradAvg(samples.getDimensions());
    for (size_t i = 0; i < samples.getSize(); ++i) {
      gradAvg.add(samples.getValues()[i]);
    }
    gradAvg.mult(1.0 / (samples.getSize()));

    for (size_t d = 0; d < samples.getDimensions(); ++d) {
      sgpp::base::DataVector col = gradAvg;
      col.mult(gradAvg[d]);
      mat.setColumn(d, col);
    }
  }

  sgpp::base::DataMatrix eigenVectors(dimensions, dimensions);
  sgpp::base::DataVector eigenValues(dimensions);
  EigenHelper::svd(EigenHelper::toEigen(mat), eigenVectors, eigenValues);
  output.eigenValues = eigenValues;
  output.eigenVectors = eigenVectors;

  double sum = eigenValues.sum();
  double partialSum = 0;
  double max = 0;
  size_t cutoff = redConfig.maxReducedDimension;
  for (size_t d = 0; d < redConfig.maxReducedDimension - 1; ++d) {
    partialSum += eigenValues[d];
    if (partialSum / sum > 0.95) {
      cutoff = d + 1;
      break;
    }
  }
  eigenVectors.transpose();
  eigenVectors.resizeRows(cutoff);
  eigenVectors.transpose();
  return eigenVectors;
}

GridReductionResult activeSubspaceReductionStep(PointSample<double>& asSample,
                                                const DataMatrix& basis,
                                                DimReduction::RegressionConfig& config) {
  SGridSample s;
  double error;
  createRegressionSurrogate(asSample, config, s, error);
  std::shared_ptr<VectorFunction> projection =
      std::make_shared<InputProjectionFunction>(basis, basis.getNcols());
  std::shared_ptr<ScalarFunction> reducedFunc = std::make_shared<InterpolantScalarFunction>(s);
  return GridReductionResult{projection, reducedFunc, s};
}

DataMatrix createBestBasis(PointSample<double>& originalSamples, PointSample<double>& samples,
    PointSample<DataVector>& gradients,
    std::vector<std::shared_ptr<ScalarFunction>>& funcs, DimReduction::BasisConfig& bConfig,
                           DimReduction::ReductionConfig& redConfig,
                           DimReduction::ExaminationConfig& examConfig,
                           DimReduction::Output& output) {
  double bestError = std::numeric_limits<double>::infinity();
  DataMatrix bestBasis;
  for (size_t i = 0; i < bConfig.basisIterations; ++i) {
    auto mat = createBasis(gradients, bConfig, redConfig, output);

    std::tuple<PointSample<double>, PointSample<double>> split = createActiveSubspaceSample(
        originalSamples, samples, mat, bConfig.evaluationConfig.sampleCount);
    auto asSample = std::get<0>(split);
    auto errorSample = std::get<1>(split);

    sgpp::base::GridReductionResult result =
        activeSubspaceReductionStep(asSample, mat, bConfig.evaluationConfig);

    auto tempFuncs = funcs;
    tempFuncs.push_back(result.replacementFunction);
    std::shared_ptr<ScalarFunction> tempSum = std::make_shared<SumScalarFunction>(tempFuncs);
    double err = examConfig.error.calculateError(tempSum, errorSample);
    std::cout << "Trying basis with dimensions: " << mat.getNcols()
              << ", " + examConfig.error.getName() + ": " << err << ".\n";
    if (err < bestError) {
      bestError = err;
      bestBasis = mat;
    }
  }
  std::cout << "Found best basis dimensions: " << bestBasis.getNcols() << " with error "
            << bestError << ".\n";
  return bestBasis;
}

void updateSample(PointSample<double>& sample, std::shared_ptr<ScalarFunction> func) {
  for (int j = 0; j < sample.getSize(); j++) {
    sample.getValues()[j] -= func->eval(sample.getKeys()[j]);
  }
}

AsReductionResult DimReduction::reduceAS(PointSample<double>& sample, ReductionConfig& redConfig,
                                         RegressionConfig& regConfig, GradientConfig& gradConfig,
                                         BasisConfig& basisConfig, ExaminationConfig& examConfig,
                                         std::vector<Output>& output) {
  auto gradientFunc = gradConfig.gradientFunction;
  std::vector<GridReductionResult> results;
  std::vector<std::shared_ptr<ScalarFunction>> funcs;
  PointSample<double> current = sample;
  std::vector<int> indices(sample.getSize());
  for (int j = 0; j < sample.getSize(); j++) indices[j] = j;
  std::random_shuffle(indices.begin(), indices.end());

  double currentError = std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < redConfig.maxIterations; i++) {
    std::cout << std::endl;
    std::cout << "Starting iteration: " << i << std::endl;

    output.push_back({});
    PointSample<DataVector> gradients = createGradientSample(current, gradConfig);

    DataMatrix basis =
        createBestBasis(sample, current, gradients, funcs, basisConfig, redConfig, examConfig, output[i]);

    std::tuple<PointSample<double>, PointSample<double>> split =
        createActiveSubspaceSample(sample, current, basis, regConfig.sampleCount);
    auto asSample = std::get<0>(split);
    auto errorSample = std::get<1>(split);

    sgpp::base::GridReductionResult result =
        activeSubspaceReductionStep(asSample, basis, regConfig);

    basisConfig.type = BasisConfig::RANDOM;

    auto tempFuncs = funcs;
    tempFuncs.push_back(result.replacementFunction);
    std::shared_ptr<ScalarFunction> tempSum = std::make_shared<SumScalarFunction>(tempFuncs);
    double err = examConfig.error.calculateError(tempSum, errorSample);
    if (!examConfig.discardWorseIterations || err < currentError) {
      {
        currentError = err;
        results.push_back(result);
        funcs.push_back(result.replacementFunction);
        updateSample(current, result.replacementFunction);
      }

      if (gradConfig.type == GradientConfig::GRADIENT_FUNCTION) {
        auto grad = finiteDifferencesFunction(tempSum, gradConfig.h);
        std::vector<std::shared_ptr<VectorFunction>> gradFuncs = {gradientFunc,
                                                                  grad};
        std::vector<bool> gradSigns = {true, false};
        std::shared_ptr<VectorFunction> gradientSum =
            std::make_shared<SumVectorFunction>(gradFuncs, gradSigns);
        gradientFunc = gradientSum;
      }

      std::cout << "-----"
                << "\n";
      std::cout << "New " + examConfig.error.getName() + " after iteration: " << currentError
                << ".\n";
    } else {
      std::cout << "Found worse " + examConfig.error.getName() + ": " << err << ".\n";
    }
  }
  std::shared_ptr<ScalarFunction> sum = std::make_shared<SumScalarFunction>(funcs);
  return AsReductionResult{sum, results};
}

ReductionResult::ReductionResult(std::shared_ptr<ScalarFunction> replacement)
    : replacementFunction(replacement) {}

double DimReduction::L2::calculateError(std::shared_ptr<ScalarFunction>& func,
                                        PointSample<double>& errorSample) {
  double res = 0;
  for (size_t i = 0; i < errorSample.getSize(); i++) {
    double val = errorSample.getValues()[i];
    double newVal = func->eval(errorSample.getKeys()[i]);
    res += pow(val - newVal, 2);
  }

  return sqrt(res / static_cast<double>(errorSample.getSize()));
}

double DimReduction::MSE::calculateError(std::shared_ptr<ScalarFunction>& func,
                                         PointSample<double>& errorSample) {
  return std::pow(L2_ERROR.calculateError(func, errorSample), 2.0);
}

double DimReduction::NRMSE::calculateError(std::shared_ptr<ScalarFunction>& func,
                                           PointSample<double>& errorSample) {
  double denom = *std::max_element(errorSample.getValues().begin(), errorSample.getValues().end()) -
                 *std::min_element(errorSample.getValues().begin(), errorSample.getValues().end());
  return L2_ERROR.calculateError(func, errorSample) / denom;
}

GridReductionResult::GridReductionResult(std::shared_ptr<VectorFunction>& t,
                                         std::shared_ptr<ScalarFunction>& sampleFunction,
                                         SGridSample& sample)
    : ReductionResult(std::make_shared<ChainScalarFunction>(
          std::vector<std::shared_ptr<VectorFunction>>(1, t), sampleFunction)),
      reducedSample(sample),
      reducedFunction(sampleFunction),
      transformation(t) {}

AsReductionResult::AsReductionResult(std::shared_ptr<ScalarFunction> replacement,
                                     std::vector<GridReductionResult>& results)
    : ReductionResult(replacement), reductions(results) {}

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
    if (sum != 0.0) {
      col.mult(1 / sum);
    }
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

    if (std::isnan(value[d])) {
      throw std::logic_error("Nan");
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


}  // namespace base
}  // namespace sgpp
