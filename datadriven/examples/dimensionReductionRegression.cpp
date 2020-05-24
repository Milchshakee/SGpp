// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_learnerRegressionTest_cpp Regression Learner
 * This example demonstrates sparse grid regression learning.
 */

#include <exception>
#include <limits>
#include <ostream>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegressionLearner.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <sgpp/datadriven/tools/dimension/DimReduction.hpp>
#include <sgpp/globaldef.hpp>
#include "sgpp/base/function/scalar/WrapperScalarFunction.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/tools/Sample.hpp"
#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/base/function/scalar/ChainScalarFunction.hpp>
#include <sgpp/base/function/vector/BoundingBoxFunction.hpp>

/**
 * @brief getLearner
 * @param dimension is the number of dimensions
 * @param regularizationConfig
 * @return a sparse grid regression learner
 */
sgpp::datadriven::RegressionLearner getLearner(
    size_t dimension, sgpp::datadriven::RegularizationConfiguration regularizationConfig) {
  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = dimension;
  gridConfig.level_ = 2;
  //  gridConfig.type_ = sgpp::base::GridType::ModLinear;

  gridConfig.type_ = sgpp::base::GridType::LinearBoundary;
  gridConfig.maxDegree_ = 3;

  auto adaptivityConfig = sgpp::base::AdaptivityConfiguration();
  adaptivityConfig.numRefinementPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-8;
  solverConfig.threshold_ = 1e-5;

  return sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                             solverConfig, regularizationConfig);
}

/**
 * @brief showRegularizationConfiguration
 * @param regularizationConfig
 * @return type of the regularization method as string
 */
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

/**
 * @brief gridSearch performs a hyper-parameter grid search over configs using a
 * holdout validation set.
 * @param configs are the regularization configs that will be tried
 * @param dimension is the number of dimensions
 * @param xTrain are the training predictors
 * @param yTrain is the training target
 * @param xValidation are the validation predictors
 * @param yValidation is the validation target
 * @return best found regularization configuration
 */
sgpp::datadriven::RegularizationConfiguration gridSearch(
    std::vector<sgpp::datadriven::RegularizationConfiguration> configs, size_t dimension,
    sgpp::base::DataMatrix& xTrain, sgpp::base::DataVector& yTrain,
    sgpp::base::DataMatrix& xValidation, sgpp::base::DataVector& yValidation) {
  double bestMSE = std::numeric_limits<double>::max();
  sgpp::datadriven::RegularizationConfiguration bestConfig;
  for (const auto& config : configs) {
    // Step 1: Create a learner
    auto learner = getLearner(dimension, config);
    // Step 2: Train it with the hyperparameter
    learner.train(xTrain, yTrain);
    // Step 3: Evaluate accuracy
    const double curMSE = learner.getMSE(xValidation, yValidation);
    std::cout << "Tested parameters are\n" << showRegularizationConfiguration(config) << ".\n";
    if (curMSE < bestMSE) {
      std::cout << "Better! RMSE is now " << std::sqrt(curMSE) << std::endl;
      bestConfig = config;
      bestMSE = curMSE;
    } else {
      std::cout << "Worse!  RMSE is now " << std::sqrt(curMSE) << std::endl;
    }
  }
  std::cout << "gridSearch finished with parameters " << showRegularizationConfiguration(bestConfig)
            << std::endl;
  return bestConfig;
}

/**
 * @brief getConfigs
 * @return some regularization configurations for seven lambdas between 1 and 0.000001
 * and for exponent bases 1.0, 0.5, 0.25, 0.125
 */
std::vector<sgpp::datadriven::RegularizationConfiguration> getConfigs() {
  decltype(getConfigs()) result;
  std::vector<double> lambdas = {1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};

  std::vector<double> exponentBases = {1.0, 0.5, 0.25, 0.125};
  for (const auto lambda : lambdas) {
    // Identity
    const auto regularizationType = sgpp::datadriven::RegularizationType::Identity;
    auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
    regularizationConfig.type_ = regularizationType;
    regularizationConfig.lambda_ = lambda;
    regularizationConfig.exponentBase_ = 0.25;
    result.push_back(regularizationConfig);
        {
          // Laplace
          const auto regularizationType = sgpp::datadriven::RegularizationType::Laplace;
          auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
          regularizationConfig.type_ = regularizationType;
          regularizationConfig.lambda_ = lambda;
          regularizationConfig.exponentBase_ = 0.25;
          result.push_back(regularizationConfig);
        }
    // Diagonal
    for (const auto exponentBase : exponentBases) {
      const auto regularizationType = sgpp::datadriven::RegularizationType::Diagonal;
      auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
      regularizationConfig.type_ = regularizationType;
      regularizationConfig.lambda_ = lambda;
      regularizationConfig.exponentBase_ = exponentBase;
      result.push_back(regularizationConfig);
    }
  }

  return result;
}

double f(const sgpp::base::DataVector& v) { return 3 * v[0] + 1 * v[1]; }

double fsin = std::sin(0.15 * M_PI);
double fcos = std::cos(0.15 * M_PI);

double g(const sgpp::base::DataVector& v) {
  double x = fcos * v[0] - fsin * v[1];
  double y = fsin * v[0] + fcos * v[1];
  return std::max(1 - std::abs(5 * x - 2.5), 0.0);
}

double ebolaFunc(const sgpp::base::DataVector& x) {
  return (x[0] + ((x[1])) + (x[2]));
}


/**
 * @brief main is an example for the RegressionLearner. It performs a grid search for the best
 * hyper-parameter for the Friedman3 dataset using the diagonal Tikhonov regularization method.
 */
int main(int argc, char** argv) {
  std::shared_ptr<sgpp::base::ScalarFunction> func =
      std::make_shared<sgpp::base::WrapperScalarFunction>(3, ebolaFunc);
  sgpp::base::BoundingBox liberiaBb({{0.1, 0.4},        // beta_1
                                     {0.1, 0.4},        // beta_2
                                     {0.05, 0.2},  // gamma_1
    });   // psi
  std::vector<sgpp::base::DistributionType> types(3, sgpp::base::DistributionType::Uniform);
  std::shared_ptr<sgpp::base::VectorFunction> bbTrans =
      std::make_shared<sgpp::base::BoundingBoxFunction>(
          sgpp::base::BoundingBoxFunction::Type::FROM_UNIT_BB, liberiaBb);
  auto vas = {bbTrans};
  std::shared_ptr<sgpp::base::ScalarFunction> unitFunc =
      std::make_shared<sgpp::base::ChainScalarFunction>(vas, func);

  std::shared_ptr<sgpp::base::DistributionUniform> u =
      std::make_shared<sgpp::base::DistributionUniform>();
  sgpp::base::DistributionsVector v(3, u);
  auto dist = sgpp::base::DistributionSample(1000, v);
  sgpp::base::ActiveSubspaceInfo i = sgpp::base::DimReduction::activeSubspaceMC(*unitFunc, dist);

  int dims = 3;
    sgpp::base::DataMatrix zonotope(dims, dims);
  for (size_t d = 0; d < dims; ++d) {
    sgpp::base::DataVector unit(dims, 0.0);
    unit[d] = 1;
    zonotope.setColumn(d, unit);
  }
  i.eigenVectors = zonotope;

  auto dist2 = sgpp::base::DistributionSample(1000, v);
  sgpp::base::PointSample<double> sample =
      sgpp::base::DimReduction::createActiveSubspaceSample(sgpp::base::SampleHelper::sampleScalarFunction(dist2, *unitFunc), i.eigenVectors, 3);
  sgpp::datadriven::Dataset data = sgpp::base::SampleHelper::fromPointSample(sample);

  
  auto dist3 = sgpp::base::DistributionSample(1000, v);
  sgpp::base::PointSample<double> valSample = sgpp::base::SampleHelper::sampleScalarFunction(dist3, *unitFunc);
  sgpp::datadriven::Dataset valData = sgpp::base::SampleHelper::fromPointSample(valSample);


  auto dataTrain = data;
  auto xTrain = dataTrain.getData();
  auto yTrain = dataTrain.getTargets();
  const auto dimensions = dataTrain.getDimension();

  auto dataValidation = valData;
  auto xValidation = dataValidation.getData();
  auto yValidation = dataValidation.getTargets();

  const auto configs = getConfigs();
  const auto bestConfig = gridSearch(configs, dimensions, xTrain, yTrain, xValidation, yValidation);

  auto dataTest = valData;
  auto xTest = dataTest.getData();
  auto yTest = dataTest.getTargets();

  auto learner = getLearner(dimensions, bestConfig);
  learner.train(xTrain, yTrain);
  const auto MSETest = learner.getMSE(xTest, yTest);
  std::cout << "Best config got a testing MSE of " << MSETest << "!" << std::endl;
}
