#pragma once

#include <sgpp/base/function/vector/VectorFunction.hpp>

namespace sgpp
{
namespace base
{
class SumVectorFunction : public VectorFunction
{

public:
  SumVectorFunction(std::vector<std::shared_ptr<VectorFunction>>& functions,
                    std::vector<bool> signs = {});

  void eval(const DataVector& x, DataVector& out) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;

private:
  std::vector<std::shared_ptr<VectorFunction>> functions;
 std::vector<bool> signs;
};
}
}
