#include "ANOVAReducer.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/hash/OperationHierarchisation.hpp"
#include "sgpp/datadriven/operation/hash/OperationPiecewiseConstantRegression/OperationPiecewiseConstantRegression.hpp"
#include <iostream>

ANOVAReducer::ANOVAReducer(size_t anovaOrder) : FunctionReducer<ANOVAInformation>(nullptr){}


void ANOVAReducer::evaluateFunction(sgpp::optimization::ScalarFunction& input,
                                    AnovaInformation& out)
{
  
}

std::unique_ptr<sgpp::optimization::ScalarFunction> ANOVAReducer::reduce(
    sgpp::optimization::ScalarFunction& input, size_t n, const AnovaInformation& info) {
  size_t dim = input.getNumberOfParameters();
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createAnovaBoundaryGrid(dim));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

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
  return nullptr;
}