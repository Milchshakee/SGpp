// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ASREDUCER_HPP
#define ASREDUCER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/tools/EigenHelper.hpp>
#include <sgpp/base/tools/dimension/DimReduction.hpp>


namespace sgpp {
namespace base {

struct AsInfo {
  sgpp::base::DataMatrix eigenVectors;
  sgpp::base::DataVector eigenValues;
};

}  // namespace base
}  // namespace sgpp

#endif
