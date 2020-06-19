#pragma once
#include <sgpp/base/function/scalar/ScalarFunction.hpp>

namespace sgpp
{
namespace base
{
class SumScalarFunction : public ScalarFunction
{

public:
  SumScalarFunction(std::vector<std::shared_ptr<ScalarFunction>>& functions, std::vector<bool> signs = {});

  double eval(const DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;

private:
  std::vector<std::shared_ptr<ScalarFunction>> functions;
 std::vector<bool> signs;
};
}
}
