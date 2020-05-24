#include <sgpp/base/tools/AnovaAnalyser.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/operation/hash/OperationEvalAnovaPrewaveletBoundary.hpp>
#include <sgpp/datadriven/tools/dimension/DimReduction.hpp>


sgpp::base::AnovaAnalyser::TransformationFunction::~TransformationFunction() {
}

void sgpp::base::AnovaAnalyser::TransformationFunction::
eval(const DataVector& x, DataVector& value) {
}

void sgpp::base::AnovaAnalyser::TransformationFunction::clone(
  std::unique_ptr<VectorFunction>& clone) const {
}

sgpp::base::SGridSample sgpp::base::AnovaAnalyser::createReducedAnovaSample(
    sgpp::base::SGridSample& sample, size_t reducedDims) {
  std::shared_ptr<sgpp::base::Grid> grid(
      sgpp::base::Grid::createAnovaPrewaveletBoundaryGrid(reducedDims));
  AnovaTypes::level_t level =
      AnovaTypes::fromNormalLevel(const_cast<Grid&>(sample.getGrid()).getStorage().getMaxLevel());
  grid->getGenerator().regular(level);
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

double calcVariance(sgpp::base::SGridSample& sample, sgpp::base::DistributionSample& dist,
                    sgpp::base::AnovaTypes::AnovaComponent& comp)
{
  std::unique_ptr<sgpp::base::OperationEval> evalOp = std::make_unique<sgpp::base::OperationEvalAnovaPrewaveletBoundary>(
      const_cast<sgpp::base::Grid&>(sample.getGrid()).getStorage(), comp);
  sgpp::base::InterpolantScalarFunction f(
      sample, evalOp);
  return 0;
}

void iterateComponentsRec(std::map<sgpp::base::AnovaTypes::AnovaComponent, double>& variances,
                          sgpp::base::SGridSample& sample, sgpp::base::DistributionSample& dist,
                          sgpp::base::AnovaTypes::AnovaComponent& comp, size_t dim) {
  if (dim == sample.getDimensions()) {
    double v = calcVariance(sample, dist, comp);
    variances.emplace(comp, v);
    return;
  }

  sgpp::base::AnovaTypes::AnovaComponent copy = comp;
  copy[dim] = false;
  iterateComponentsRec(variances, sample, dist, copy, dim + 1);
  copy[dim] = true;
  iterateComponentsRec(variances, sample, dist, copy, dim + 1);
}

sgpp::base::AnovaAnalyser::Info sgpp::base::AnovaAnalyser::evaluateFunction(
    sgpp::base::SGridSample& sample, DistributionSample& dist) {
  AnovaTypes::AnovaComponent comp(sample.getDimensions());
  std::map<sgpp::base::AnovaTypes::AnovaComponent, double> variances;
  iterateComponentsRec(variances, sample, dist, comp, 0);
    }

sgpp::base::AnovaAnalyser::TransformationFunction::TransformationFunction(
    size_t d, const std::vector<size_t>& dimension_order, size_t reduced_dims)
    : VectorFunction(d, reduced_dims), dimensionOrder(dimension_order) {}

size_t sgpp::base::AnovaAnalyser::getReducedDimensions(Info& i, double minVariance) {
  double covered = 1;
  for (int j = i.dimensionOrder.size() - 1; j >= 0; --j) {
    covered -= i.totalEffectIndices[j];
    if (covered < minVariance) {
      return j + 1;
    }
  }
  return 1;
}
