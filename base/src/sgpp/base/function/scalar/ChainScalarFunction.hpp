#pragma once

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>

namespace sgpp {
namespace base {

class ChainScalarFunction : public ScalarFunction
{

public:
  ChainScalarFunction(const std::vector<std::shared_ptr<VectorFunction>>& chain,
                      const std::shared_ptr < ScalarFunction>& func);
  ~ChainScalarFunction() override;
  double eval(const DataVector& x) override;
  void clone(std::unique_ptr<ScalarFunction>& clone) const override;
private:
  std::vector<std::shared_ptr<VectorFunction>> chain;
 std::shared_ptr < ScalarFunction> func;
};
}  // namespace base
}  // namespace sgpp