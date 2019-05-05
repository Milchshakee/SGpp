#include "ANOVAReducer.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/operation/hash/OperationHierarchisation.hpp"

ANOVAReducer::ANOVAReducer(size_t anovaOrder) {}

std::unique_ptr<sgpp::optimization::ScalarFunction> ANOVAReducer::reduceFunction(
    sgpp::optimization::ScalarFunction& input) {
  size_t dim = input.getNumberOfParameters();
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 1));
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
  //std::unique_ptr<sgpp::base::OperationHierarchisation>(
  //    sgpp::op_factory::createOperationHierarchisation(*grid))
  //    ->doHierarchisation(alpha);

  //sgpp::base::DataVector p(dim);
  //p[0] = 0.52;
  //p[1] = 0.73;
  //std::unique_ptr<sgpp::base::OperationEval> opEval(sgpp::op_factory::createOperationEval(*grid));
  //std::cout << "u(0.52, 0.73) = " << opEval->eval(alpha, p) << std::endl;
  return nullptr;
}
