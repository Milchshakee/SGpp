#include <iostream>
#include "ANOVAReducer.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/hash/OperationHierarchisation.hpp"
#include <sgpp/base/operation/BaseOpFactory.hpp>

ANOVAReducer::ANOVAReducer(size_t gridLevel) : FunctionReducer<A>(nullptr) {}

void ANOVAReducer::evaluateFunction(sgpp::optimization::ScalarFunction& input,
                                    AnovaInformation& out)
{
  size_t dim = input.getNumberOfParameters();
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createAnovaBoundaryGrid(dim));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(gridLevel);

  sgpp::base::DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector arg(dim);

    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      arg[j] = gp.getStandardCoordinate(j);
    }
    alpha[i] = input.eval(arg);
  }

  std::unique_ptr<sgpp::base::OperationHierarchisation>(
      sgpp::op_factory::createOperationHierarchisation(*grid))
      ->doHierarchisation(alpha);
  std::unique_ptr<sgpp::base::OperationAnova> op(sgpp::op_factory::createOperationAnova(*grid));

  out.variances = std::vector<double>(input.getNumberOfParameters(), 0.0);
  for (size_t d = 0; d < input.getNumberOfParameters(); ++d) {
    std::vector<sgpp::base::OperationAnova::AnovaComponent> comps =
        op->calculateAnovaOrderVariance(alpha, d);
    for (sgpp::base::OperationAnova::AnovaComponent c : comps) {
      out.variances[d] += c.variance;
      out.components.push_back({c.order, c.variance, c.fixedDimensions});
    }
  }
}

std::unique_ptr<sgpp::optimization::ScalarFunction> ANOVAReducer::reduce(
    sgpp::optimization::ScalarFunction& input, size_t n, const AnovaInformation& info) {

  return nullptr;
}