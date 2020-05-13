#pragma once
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

namespace sgpp
{
namespace base
{
class BoundingBoxFunction : public VectorFunction
{
public:
  enum Type { TO_UNIT_BB, FROM_UNIT_BB
  };

  BoundingBoxFunction(Type type, const BoundingBox& bb);
  void eval(const DataVector& x, DataVector& value) override;
  void clone(std::unique_ptr<VectorFunction>& clone) const override;

private:
  Type type;
  BoundingBox bb;
};
}
}
