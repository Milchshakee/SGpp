// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearGridStencil.hpp>
#include <sgpp/base/grid/type/LinearL0BoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModLinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGridStencil.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/ModWaveletGrid.hpp>
#include <sgpp/base/grid/type/PeriodicGrid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/grid/type/PolyClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/PolyClenshawCurtisBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModPolyClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryCombigridGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>
#include <sgpp/base/grid/type/SquareRootGrid.hpp>
#include <sgpp/base/grid/type/WaveletBoundaryGrid.hpp>
#include <sgpp/base/grid/type/WaveletGrid.hpp>
#include <sgpp/base/grid/type/AnovaLinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/AnovaPrewaveletBoundaryGrid.hpp>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp>
#include <sgpp/globaldef.hpp>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace base {

Grid* Grid::createLinearGridStencil(size_t dim) { return new LinearGridStencil(dim); }

Grid* Grid::createModLinearGridStencil(size_t dim) { return new ModLinearGridStencil(dim); }

Grid* Grid::createLinearGrid(size_t dim) { return new LinearGrid(dim); }

Grid* Grid::createLinearStretchedGrid(size_t dim) { return new LinearStretchedGrid(dim); }

Grid* Grid::createAnovaLinearBoundaryGrid(size_t dim) {
  return new AnovaLinearBoundaryGrid(dim);
}

Grid* Grid::createAnovaLinearBoundaryGrid(size_t dim,
                                          std::vector<AnovaTypes::LevelIndexPair>& anchor) {
  return new AnovaLinearBoundaryGrid(dim, anchor);
}

  Grid* Grid::createAnovaPrewaveletBoundaryGrid(size_t dim) {
  return new AnovaPrewaveletBoundaryGrid(dim);
}

    Grid* Grid::createAnovaPrewaveletBoundaryGrid(size_t dim, std::vector<AnovaTypes::LevelIndexPair>& anchor) {
  return new AnovaPrewaveletBoundaryGrid(dim, anchor);
}

Grid* Grid::createLinearBoundaryGrid(size_t dim, level_t boundaryLevel) {
  if (boundaryLevel == 0) {
    return new LinearL0BoundaryGrid(dim);
  } else {
    return new LinearBoundaryGrid(dim, boundaryLevel);
  }
}

Grid* Grid::createLinearStretchedBoundaryGrid(size_t dim) {
  return new LinearStretchedBoundaryGrid(dim);
}

Grid* Grid::createLinearClenshawCurtisGrid(size_t dim) { return new LinearClenshawCurtisGrid(dim); }

Grid* Grid::createLinearClenshawCurtisBoundaryGrid(size_t dim, level_t boundaryLevel) {
  return new LinearClenshawCurtisBoundaryGrid(dim, boundaryLevel);
}

Grid* Grid::createModLinearClenshawCurtisGrid(size_t dim) {
  return new ModLinearClenshawCurtisGrid(dim);
}

Grid* Grid::createModLinearGrid(size_t dim) { return new ModLinearGrid(dim); }

Grid* Grid::createPolyGrid(size_t dim, size_t degree) { return new PolyGrid(dim, degree); }

Grid* Grid::createPolyBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel) {
  return new PolyBoundaryGrid(dim, degree, boundaryLevel);
}

Grid* Grid::createPolyClenshawCurtisGrid(size_t dim, size_t degree) {
  return new PolyClenshawCurtisGrid(dim, degree);
}

Grid* Grid::createModPolyClenshawCurtisGrid(size_t dim, size_t degree) {
  return new ModPolyClenshawCurtisGrid(dim, degree);
}

Grid* Grid::createPolyClenshawCurtisBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel) {
  return new PolyClenshawCurtisBoundaryGrid(dim, degree, boundaryLevel);
}

Grid* Grid::createWaveletGrid(size_t dim) { return new WaveletGrid(dim); }

Grid* Grid::createWaveletBoundaryGrid(size_t dim) { return new WaveletBoundaryGrid(dim); }

Grid* Grid::createModWaveletGrid(size_t dim) { return new ModWaveletGrid(dim); }

Grid* Grid::createBsplineGrid(size_t dim, size_t degree) { return new BsplineGrid(dim, degree); }

Grid* Grid::createBsplineBoundaryGrid(size_t dim, size_t degree) {
  return new BsplineBoundaryGrid(dim, degree);
}

Grid* Grid::createBsplineClenshawCurtisGrid(size_t dim, size_t degree) {
  return new BsplineClenshawCurtisGrid(dim, degree);
}

Grid* Grid::createModBsplineGrid(size_t dim, size_t degree) {
  return new ModBsplineGrid(dim, degree);
}

Grid* Grid::createModBsplineClenshawCurtisGrid(size_t dim, size_t degree) {
  return new ModBsplineClenshawCurtisGrid(dim, degree);
}

Grid* Grid::createFundamentalSplineGrid(size_t dim, size_t degree) {
  return new FundamentalSplineGrid(dim, degree);
}

Grid* Grid::createModFundamentalSplineGrid(size_t dim, size_t degree) {
  return new ModFundamentalSplineGrid(dim, degree);
}

Grid* Grid::createSquareRootGrid(size_t dim) { return new SquareRootGrid(dim); }

Grid* Grid::createPrewaveletGrid(size_t dim) { return new PrewaveletGrid(dim); }

Grid* Grid::createLinearTruncatedBoundaryGrid(size_t dim) {
  return new LinearTruncatedBoundaryGrid(dim);
}

Grid* Grid::createModPolyGrid(size_t dim, size_t degree) { return new ModPolyGrid(dim, degree); }

Grid* Grid::createPeriodicGrid(size_t dim) { return new PeriodicGrid(dim); }

Grid* Grid::createNakBsplineBoundaryCombigridGrid(size_t dim, size_t degree) {
  return new NakBsplineBoundaryCombigridGrid(dim, degree);
}

Grid* Grid::createGrid(RegularGridConfiguration gridConfig) {
  if (gridConfig.filename_.length() > 0) {
    std::ifstream ifs(gridConfig.filename_);
    std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    return base::Grid::unserialize(content);
  } else {
    switch (gridConfig.type_) {
      case GridType::AnovaLinearBoundary:
        return Grid::createAnovaLinearBoundaryGrid(gridConfig.dim_);
      case GridType::AnovaPrewaveletBoundary:
        return Grid::createAnovaPrewaveletBoundaryGrid(gridConfig.dim_);
      case GridType::Linear:
        return Grid::createLinearGrid(gridConfig.dim_);
      case GridType::LinearStretched:
        return Grid::createLinearStretchedGrid(gridConfig.dim_);
      case GridType::LinearL0Boundary:
        return Grid::createLinearBoundaryGrid(gridConfig.dim_);
      case GridType::LinearBoundary:
        return Grid::createLinearBoundaryGrid(gridConfig.dim_, gridConfig.boundaryLevel_);
      case GridType::LinearStretchedBoundary:
        return Grid::createLinearStretchedBoundaryGrid(gridConfig.dim_);
      case GridType::LinearTruncatedBoundary:
        return Grid::createLinearTruncatedBoundaryGrid(gridConfig.dim_);
      case GridType::ModLinear:
        return Grid::createModLinearGrid(gridConfig.dim_);
      case GridType::Poly:
        return Grid::createPolyGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::PolyBoundary:
        return Grid::createPolyBoundaryGrid(gridConfig.dim_, gridConfig.maxDegree_,
                                            gridConfig.boundaryLevel_);
      case GridType::ModPoly:
        return Grid::createModPolyGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::PolyClenshawCurtis:
        return Grid::createPolyClenshawCurtisGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::PolyClenshawCurtisBoundary:
        return Grid::createPolyClenshawCurtisBoundaryGrid(gridConfig.dim_, gridConfig.maxDegree_,
                                                          gridConfig.boundaryLevel_);
      case GridType::ModPolyClenshawCurtis:
        return Grid::createModPolyClenshawCurtisGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::ModWavelet:
        return Grid::createModWaveletGrid(gridConfig.dim_);
      case GridType::ModBspline:
        return Grid::createModBsplineGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::Prewavelet:
        return Grid::createPrewaveletGrid(gridConfig.dim_);
      case GridType::SquareRoot:
        return Grid::createSquareRootGrid(gridConfig.dim_);
      case GridType::Periodic:
        return Grid::createPeriodicGrid(gridConfig.dim_);
      case GridType::LinearClenshawCurtisBoundary:
        return Grid::createLinearClenshawCurtisBoundaryGrid(gridConfig.dim_,
                                                            gridConfig.boundaryLevel_);
      case GridType::LinearClenshawCurtis:
        return Grid::createLinearClenshawCurtisGrid(gridConfig.dim_);
      case GridType::ModLinearClenshawCurtis:
        return Grid::createModLinearClenshawCurtisGrid(gridConfig.dim_);
      case GridType::Bspline:
        return Grid::createBsplineGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::BsplineBoundary:
        return Grid::createBsplineBoundaryGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::BsplineClenshawCurtis:
        return Grid::createBsplineClenshawCurtisGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::Wavelet:
        return Grid::createWaveletGrid(gridConfig.dim_);
      case GridType::WaveletBoundary:
        return Grid::createWaveletBoundaryGrid(gridConfig.dim_);
      case GridType::FundamentalSpline:
        return Grid::createFundamentalSplineGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::ModFundamentalSpline:
        return Grid::createModFundamentalSplineGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::ModBsplineClenshawCurtis:
        return Grid::createModBsplineClenshawCurtisGrid(gridConfig.dim_, gridConfig.maxDegree_);
      case GridType::LinearStencil:
        return Grid::createLinearStretchedGrid(gridConfig.dim_);
      case GridType::ModLinearStencil:
        return Grid::createModLinearGridStencil(gridConfig.dim_);
      case GridType::NakBsplineBoundaryCombigrid:
        return Grid::createNakBsplineBoundaryCombigridGrid(gridConfig.dim_, gridConfig.maxDegree_);
      default:
        throw generation_exception("Grid::createGrid - grid type not known");
    }
  }
}

Grid* Grid::createGridOfEquivalentType(size_t numDims) {
  Grid* newGrid = nullptr;

  level_t boundaryLevel = 0;
  size_t degree = 1;
  switch (getType()) {
    case GridType::Linear:
      newGrid = Grid::createLinearGrid(numDims);
      break;
    case GridType::LinearStretched:
      newGrid = Grid::createLinearStretchedGrid(numDims);
      break;
    case GridType::LinearL0Boundary:
      newGrid = Grid::createLinearBoundaryGrid(numDims);
      break;
    case GridType::AnovaLinearBoundary:
      newGrid = Grid::createAnovaLinearBoundaryGrid(numDims);
      break;
    case GridType::AnovaPrewaveletBoundary:
      newGrid = Grid::createAnovaPrewaveletBoundaryGrid(numDims);
      break;
    case GridType::LinearBoundary:
      boundaryLevel =
          dynamic_cast<BoundaryGridGenerator*>(&this->getGenerator())->getBoundaryLevel();
      newGrid = Grid::createLinearBoundaryGrid(numDims, boundaryLevel);
      break;
    case GridType::LinearStretchedBoundary:
      newGrid = Grid::createLinearStretchedBoundaryGrid(numDims);
      break;
    case GridType::LinearTruncatedBoundary:
      newGrid = Grid::createLinearTruncatedBoundaryGrid(numDims);
      break;
    case GridType::ModLinear:
      newGrid = Grid::createModLinearGrid(numDims);
      break;
    case GridType::Poly:
      degree = dynamic_cast<PolyGrid*>(this)->getDegree();
      newGrid = Grid::createPolyGrid(numDims, degree);
      break;
    case GridType::PolyBoundary:
      degree = dynamic_cast<PolyBoundaryGrid*>(this)->getDegree();
      boundaryLevel =
          dynamic_cast<BoundaryGridGenerator*>(&this->getGenerator())->getBoundaryLevel();
      newGrid = Grid::createPolyBoundaryGrid(numDims, degree, boundaryLevel);
      break;
    case GridType::ModPoly:
      degree = dynamic_cast<ModPolyGrid*>(this)->getDegree();
      newGrid = Grid::createModPolyGrid(numDims, degree);
      break;
    case GridType::ModWavelet:
      newGrid = Grid::createModWaveletGrid(numDims);
      break;
    case GridType::ModBspline:
      degree = dynamic_cast<ModBsplineGrid*>(this)->getDegree();
      newGrid = Grid::createModBsplineGrid(numDims, degree);
      break;
    case GridType::Prewavelet:
      newGrid = Grid::createPrewaveletGrid(numDims);
      break;
    case GridType::SquareRoot:
      newGrid = Grid::createSquareRootGrid(numDims);
      break;
    case GridType::Periodic:
      newGrid = Grid::createPeriodicGrid(numDims);
      break;
    case GridType::LinearClenshawCurtisBoundary:
      boundaryLevel =
          dynamic_cast<BoundaryGridGenerator*>(&this->getGenerator())->getBoundaryLevel();
      newGrid = Grid::createLinearClenshawCurtisBoundaryGrid(numDims, boundaryLevel);
      break;
    case GridType::LinearClenshawCurtis:
      newGrid = Grid::createLinearClenshawCurtisGrid(numDims);
      break;
    case GridType::ModLinearClenshawCurtis:
      newGrid = Grid::createModLinearClenshawCurtisGrid(numDims);
      break;
    case GridType::Bspline:
      degree = dynamic_cast<BsplineGrid*>(this)->getDegree();
      newGrid = Grid::createBsplineGrid(numDims, degree);
      break;
    case GridType::BsplineBoundary:
      degree = dynamic_cast<BsplineBoundaryGrid*>(this)->getDegree();
      newGrid = Grid::createBsplineBoundaryGrid(numDims, degree);
      break;
    case GridType::BsplineClenshawCurtis:
      degree = dynamic_cast<BsplineClenshawCurtisGrid*>(this)->getDegree();
      newGrid = Grid::createBsplineClenshawCurtisGrid(numDims, degree);
      break;
    case GridType::Wavelet:
      newGrid = Grid::createWaveletGrid(numDims);
      break;
    case GridType::WaveletBoundary:
      newGrid = Grid::createWaveletBoundaryGrid(numDims);
      break;
    case GridType::FundamentalSpline:
      degree = dynamic_cast<FundamentalSplineGrid*>(this)->getDegree();
      newGrid = Grid::createFundamentalSplineGrid(numDims, degree);
      break;
    case GridType::ModFundamentalSpline:
      degree = dynamic_cast<ModFundamentalSplineGrid*>(this)->getDegree();
      newGrid = Grid::createModFundamentalSplineGrid(numDims, degree);
      break;
    case GridType::ModBsplineClenshawCurtis:
      degree = dynamic_cast<ModBsplineClenshawCurtisGrid*>(this)->getDegree();
      newGrid = Grid::createModBsplineClenshawCurtisGrid(numDims, degree);
      break;
    case GridType::LinearStencil:
      newGrid = Grid::createLinearStretchedGrid(numDims);
      break;
    case GridType::ModLinearStencil:
      newGrid = Grid::createModLinearGridStencil(numDims);
      break;
    case GridType::PolyClenshawCurtis:
      degree = dynamic_cast<PolyClenshawCurtisGrid*>(this)->getDegree();
      newGrid = Grid::createPolyClenshawCurtisGrid(numDims, degree);
      break;
    case GridType::PolyClenshawCurtisBoundary:
      degree = dynamic_cast<PolyClenshawCurtisBoundaryGrid*>(this)->getDegree();
      boundaryLevel =
          dynamic_cast<BoundaryGridGenerator*>(&this->getGenerator())->getBoundaryLevel();
      newGrid = Grid::createPolyClenshawCurtisBoundaryGrid(numDims, degree, boundaryLevel);
      break;
    case GridType::ModPolyClenshawCurtis:
      degree = dynamic_cast<ModPolyClenshawCurtisGrid*>(this)->getDegree();
      newGrid = Grid::createModPolyClenshawCurtisGrid(numDims, degree);
      break;
    case GridType::NakBsplineBoundaryCombigrid:
      degree = dynamic_cast<NakBsplineBoundaryCombigridGrid*>(this)->getDegree();
      return Grid::createNakBsplineBoundaryCombigridGrid(numDims, degree);
    default:
      throw generation_exception("Grid::clone - grid type not known");
  }
  return newGrid;
}

Grid* Grid::clone() {
  // clone grid of the same type
  Grid* newGrid = createGridOfEquivalentType(getDimension());
  newGrid->storage = this->storage;

  return newGrid;
}

GridType Grid::getZeroBoundaryType() {
  switch (getType()) {
    case GridType::Linear:
    case GridType::LinearL0Boundary:
    case GridType::LinearBoundary:
    case GridType::LinearTruncatedBoundary:
    case GridType::ModLinear:
    case GridType::SquareRoot:
    case GridType::Periodic:
    case GridType::LinearStencil:
    case GridType::ModLinearStencil:
    case GridType::AnovaLinearBoundary:
      return GridType::Linear;
    case GridType::LinearStretched:
    case GridType::LinearStretchedBoundary:
      return GridType::LinearStretched;
    case GridType::Poly:
    case GridType::PolyBoundary:
    case GridType::ModPoly:
      return GridType::Poly;
    case GridType::ModWavelet:
    case GridType::Wavelet:
    case GridType::WaveletBoundary:
    case GridType::AnovaPrewaveletBoundary:
      return GridType::Wavelet;
    case GridType::Bspline:
    case GridType::BsplineBoundary:
    case GridType::ModBspline:
      return GridType::Bspline;
    case GridType::Prewavelet:
      return GridType::Prewavelet;
    case GridType::LinearClenshawCurtis:
    case GridType::LinearClenshawCurtisBoundary:
    case GridType::ModLinearClenshawCurtis:
      return GridType::LinearClenshawCurtis;
    case GridType::FundamentalSpline:
    case GridType::ModFundamentalSpline:
      return GridType::FundamentalSpline;
    case GridType::PolyClenshawCurtis:
    case GridType::PolyClenshawCurtisBoundary:
    case GridType::ModPolyClenshawCurtis:
      return GridType::PolyClenshawCurtis;
    // no non-boundary treatment basis available for the following grids
    case GridType::BsplineClenshawCurtis:
    case GridType::ModBsplineClenshawCurtis:
    case GridType::NakBsplineBoundaryCombigrid:
    default:
      throw generation_exception("Grid::getZeroBoundaryType - no conversion known");
  }
}

std::string Grid::getTypeAsString() { return typeVerboseMap()[getType()]; }

Grid* Grid::unserialize(const std::string& istr) {
  std::istringstream istream;
  istream.str(istr);

  return Grid::unserialize(istream);
}

Grid* Grid::unserialize(std::istream& istr) {
  std::string gridtype;
  istr >> gridtype;

  if (typeMap().count(gridtype) > 0) {
    // it calls a function pointer out of a map.
    return typeMap()[gridtype](istr);
  } else {
    // compose error message and throw factory_exception
    std::string errMsg = "factory_exception unserialize: unknown gridtype.\n";
    errMsg += "Got: '" + gridtype + "'; possible choices: ";
    factoryMap::iterator it;

    for (it = typeMap().begin(); it != typeMap().end(); it++) {
      errMsg += "'" + (*it).first + "' ";
    }

    throw factory_exception(errMsg.c_str());
  }

  return nullptr;
}

std::map<std::string, Grid::Factory>& Grid::typeMap() {
  // This is only executed once!
  static factoryMap* tMap = new factoryMap();

  if (tMap->size() == 0) {
/*
 * Insert factories here. This methods may NOT read the grid type.
 * This map takes a string, function pointer pair.
 */
#ifdef _WIN32
    tMap->insert(std::pair<std::string, Grid::Factory>("NULL", Grid::nullFactory));
    tMap->insert(std::pair<std::string, Grid::Factory>("linear", LinearGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("linearStretched", LinearStretchedGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("linearL0Boundary",
                                                       LinearL0BoundaryGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("linearstencil", LinearGridStencil::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modlinearstencil",
                                                       ModLinearGridStencil::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("linearBoundary", LinearBoundaryGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("anovaLinearBoundary", AnovaLinearBoundaryGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("anovaPrewaveletBoundary", AnovaPrewaveletBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("linearStretchedBoundary",
                                                       LinearStretchedBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("linearClenshawCurtis",
                                                       LinearClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>(
        "linearClenshawCurtisBoundary", LinearClenshawCurtisBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modLinearClenshawCurtis",
                                                       ModLinearClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modlinear", ModLinearGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("poly", PolyGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("polyBoundary", PolyBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modpoly", ModPolyGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("polyClenshawCurtis",
                                                       PolyClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>(
        "polyClenshawCurtisBoundary", PolyClenshawCurtisBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modPolyClenshawCurtis",
                                                       ModPolyClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("wavelet", WaveletGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("waveletBoundary", WaveletBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modWavelet", ModWaveletGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("bspline", BsplineGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("bsplineBoundary", BsplineBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("bsplineClenshawCurtis",
                                                       BsplineClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modBspline", ModBsplineGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("fundamentalSpline",
                                                       FundamentalSplineGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modFundamentalSpline",
                                                       ModFundamentalSplineGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modBsplineClenshawCurtis",
                                                       ModBsplineClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("prewavelet", PrewaveletGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("periodic", PeriodicGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("linearTruncatedBoundary",
                                                       LinearTruncatedBoundaryGrid::unserialize));
#else
    tMap->insert(std::make_pair("NULL", Grid::nullFactory));
    tMap->insert(std::make_pair("linear", LinearGrid::unserialize));
    tMap->insert(std::make_pair("linearStretched", LinearStretchedGrid::unserialize));
    tMap->insert(std::make_pair("linearL0Boundary", LinearL0BoundaryGrid::unserialize));
    tMap->insert(std::make_pair("linearstencil", LinearGridStencil::unserialize));
    tMap->insert(std::make_pair("modlinearstencil", ModLinearGridStencil::unserialize));
    tMap->insert(std::make_pair("linearBoundary", LinearBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("anovaLinearBoundary", AnovaLinearBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("anovaPrewaveletBoundary", AnovaPrewaveletBoundaryGrid::unserialize));
    tMap->insert(
        std::make_pair("linearStretchedBoundary", LinearStretchedBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("linearClenshawCurtis", LinearClenshawCurtisGrid::unserialize));
    tMap->insert(std::make_pair("modLinearClenshawCurtis", LinearClenshawCurtisGrid::unserialize));
    tMap->insert(std::make_pair("linearClenshawCurtisBoundary",
                                LinearClenshawCurtisBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("modlinear", ModLinearGrid::unserialize));
    tMap->insert(std::make_pair("poly", PolyGrid::unserialize));
    tMap->insert(std::make_pair("polyBoundary", PolyBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("modpoly", ModPolyGrid::unserialize));
    tMap->insert(std::make_pair("polyClenshawCurtis", PolyClenshawCurtisGrid::unserialize));
    tMap->insert(
        std::make_pair("polyClenshawCurtisBoundary", PolyClenshawCurtisBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("modPolyClenshawCurtis", ModPolyClenshawCurtisGrid::unserialize));
    tMap->insert(std::make_pair("wavelet", WaveletGrid::unserialize));
    tMap->insert(std::make_pair("waveletBoundary", WaveletBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("modWavelet", ModWaveletGrid::unserialize));
    tMap->insert(std::make_pair("bspline", BsplineGrid::unserialize));
    tMap->insert(std::make_pair("bsplineBoundary", BsplineBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("bsplineClenshawCurtis", BsplineClenshawCurtisGrid::unserialize));
    tMap->insert(std::make_pair("modBspline", ModBsplineGrid::unserialize));
    tMap->insert(std::make_pair("fundamentalSpline", FundamentalSplineGrid::unserialize));
    tMap->insert(std::make_pair("modFundamentalSpline", ModFundamentalSplineGrid::unserialize));
    tMap->insert(
        std::make_pair("modBsplineClenshawCurtis", ModBsplineClenshawCurtisGrid::unserialize));
    tMap->insert(std::make_pair("prewavelet", PrewaveletGrid::unserialize));
    tMap->insert(std::make_pair("periodic", PeriodicGrid::unserialize));
    tMap->insert(
        std::make_pair("linearTruncatedBoundary", LinearTruncatedBoundaryGrid::unserialize));
#endif
  }

  return *tMap;
}

std::map<sgpp::base::GridType, std::string>& Grid::typeVerboseMap() {
  // This is only executed once!
  static gridTypeVerboseMap* verboseMap = new gridTypeVerboseMap();

  if (verboseMap->size() == 0) {
/*
 * Insert strings here.
 */
#ifdef _WIN32
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::Linear, "linear"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::LinearStretched, "linearStretched"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::LinearL0Boundary,
                                                                    "linearL0Boundary"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::LinearStencil, "linearstencil"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::ModLinearStencil,
                                                                    "modlinearstencil"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::LinearBoundary, "linearBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::AnovaLinearBoundary,
                                                                    "anovaLinearBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::AnovaPrewaveletBoundary,
                                                                    "anovaPrewaveletBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::LinearStretchedBoundary, "linearStretchedBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::LinearClenshawCurtis,
                                                                    "linearClenshawCurtis"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::LinearClenshawCurtisBoundary, "linearClenshawCurtisBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::ModLinearClenshawCurtis, "modLinearClenshawCurtis"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::ModLinear, "modlinear"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::Poly, "poly"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::PolyBoundary, "polyBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::PolyClenshawCurtis,
                                                                    "polyClenshawCurtis"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::PolyClenshawCurtisBoundary, "polyClenshawCurtisBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::ModPolyClenshawCurtis,
                                                                    "modPolyClenshawCurtis"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::ModPoly, "modpoly"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::Wavelet, "wavelet"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::WaveletBoundary, "waveletBoundary"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::ModWavelet, "modWavelet"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::Bspline, "bspline"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::BsplineBoundary, "bsplineBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::BsplineClenshawCurtis,
                                                                    "bsplineClenshawCurtis"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::ModBspline, "modBspline"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::FundamentalSpline,
                                                                    "fundamentalSpline"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::ModFundamentalSpline,
                                                                    "modFundamentalSpline"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::ModBsplineClenshawCurtis, "modBsplineClenshawCurtis"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::Prewavelet, "prewavelet"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::Periodic, "periodic"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::LinearTruncatedBoundary, "linearTruncatedBoundary"));
#else
    verboseMap->insert(std::make_pair(GridType::Linear, "linear"));
    verboseMap->insert(std::make_pair(GridType::LinearStretched, "linearStretched"));
    verboseMap->insert(std::make_pair(GridType::LinearL0Boundary, "linearL0Boundary"));
    verboseMap->insert(std::make_pair(GridType::LinearStencil, "linearstencil"));
    verboseMap->insert(std::make_pair(GridType::ModLinearStencil, "modlinearstencil"));
    verboseMap->insert(std::make_pair(GridType::LinearBoundary, "linearBoundary"));
    verboseMap->insert(std::make_pair(GridType::AnovaLinearBoundary, "anovaLinearBoundary"));
    verboseMap->insert(
        std::make_pair(GridType::AnovaPrewaveletBoundary, "anovaPrewaveletBoundary"));
    verboseMap->insert(
        std::make_pair(GridType::LinearStretchedBoundary, "linearStretchedBoundary"));
    verboseMap->insert(std::make_pair(GridType::LinearClenshawCurtis, "linearClenshawCurtis"));
    verboseMap->insert(
        std::make_pair(GridType::LinearClenshawCurtisBoundary, "linearClenshawCurtisBoundary"));
    verboseMap->insert(
        std::make_pair(GridType::ModLinearClenshawCurtis, "modLinearClenshawCurtis"));
    verboseMap->insert(std::make_pair(GridType::ModLinear, "modlinear"));
    verboseMap->insert(std::make_pair(GridType::Poly, "poly"));
    verboseMap->insert(std::make_pair(GridType::PolyBoundary, "polyBoundary"));
    verboseMap->insert(std::make_pair(GridType::ModPoly, "modpoly"));
    verboseMap->insert(std::make_pair(GridType::Wavelet, "wavelet"));
    verboseMap->insert(std::make_pair(GridType::WaveletBoundary, "waveletBoundary"));
    verboseMap->insert(std::make_pair(GridType::ModWavelet, "modWavelet"));
    verboseMap->insert(std::make_pair(GridType::Bspline, "bspline"));
    verboseMap->insert(std::make_pair(GridType::BsplineBoundary, "bsplineBoundary"));
    verboseMap->insert(std::make_pair(GridType::BsplineClenshawCurtis, "bsplineClenshawCurtis"));
    verboseMap->insert(std::make_pair(GridType::ModBspline, "modBspline"));
    verboseMap->insert(std::make_pair(GridType::FundamentalSpline, "fundamentalSpline"));
    verboseMap->insert(std::make_pair(GridType::ModFundamentalSpline, "modFundamentalSpline"));
    verboseMap->insert(
        std::make_pair(GridType::ModBsplineClenshawCurtis, "modBsplineClenshawCurtis"));
    verboseMap->insert(std::make_pair(GridType::Prewavelet, "prewavelet"));
    verboseMap->insert(std::make_pair(GridType::Periodic, "periodic"));
    verboseMap->insert(
        std::make_pair(GridType::LinearTruncatedBoundary, "linearTruncatedBoundary"));
    verboseMap->insert(std::make_pair(GridType::PolyClenshawCurtis, "polyClenshawCurtis"));
    verboseMap->insert(
        std::make_pair(GridType::PolyClenshawCurtisBoundary, "polyClenshawCurtisBoundary"));
    verboseMap->insert(std::make_pair(GridType::ModPolyClenshawCurtis, "modPolyClenshawCurtis"));
    verboseMap->insert(
        std::make_pair(GridType::ModLinearClenshawCurtis, "modLinearClenshawCurtis"));
    verboseMap->insert(std::make_pair(GridType::LinearClenshawCurtis, "linearClenshawCurtis"));
#endif
  }

  return *verboseMap;
}

/**
 * Factory for everything we don't know.
 */
Grid* Grid::nullFactory(std::istream&) {
  throw factory_exception("factory_exeception unserialize: unsupported gridtype");
  return nullptr;
}

Grid::Grid(std::istream& istr) : storage(istr) {}

Grid::Grid(size_t dim) : storage(dim) {}

Grid::Grid(BoundingBox& boundingBox) : storage(boundingBox) {}

Grid::Grid(Stretching& stretching) : storage(stretching) {}

Grid::~Grid() {}

GridStorage& Grid::getStorage() { return storage; }

BoundingBox& Grid::getBoundingBox() { return *storage.getBoundingBox(); }

Stretching& Grid::getStretching() {
  auto* stretching = storage.getStretching();

  if (stretching != nullptr) {
    return *stretching;
  } else {
    throw generation_exception("Grid does not use stretching.");
  }
}

void Grid::setBoundingBox(BoundingBox& boundingBox) { storage.setBoundingBox(boundingBox); }

void Grid::setStretching(Stretching& stretching) { storage.setStretching(stretching); }

void Grid::serialize(std::string& ostr, int version) {
  std::ostringstream ostream;
  this->serialize(ostream, version);

  ostr = ostream.str();
}

std::string Grid::serialize(int version) {
  std::ostringstream ostream;
  this->serialize(ostream, version);

  return ostream.str();
}

void Grid::serialize(std::ostream& ostr, int version) {
  ostr << typeVerboseMap()[this->getType()] << std::endl;
  storage.serialize(ostr, version);
}

void Grid::refine(DataVector& vector, int numOfPoints) {
  SurplusRefinementFunctor functor(vector, numOfPoints);
  getGenerator().refine(functor);
}

void Grid::insertPoint(size_t dim, unsigned int levels[], unsigned int indices[], bool isLeaf) {
  // create HashGridPoint object for the point
  GridPoint pointIndex(dim);

  for (unsigned int i = 0; i < dim - 1; i++) {
    pointIndex.push(i, levels[i], indices[i]);
  }

  // insert last level/index and hash
  pointIndex.set(dim - 1, levels[dim - 1], indices[dim - 1], isLeaf);
  // insert point to the GridStorage
  storage.insert(pointIndex);
}

size_t Grid::getDimension() const { return storage.getDimension(); }

size_t Grid::getSize() const { return storage.getSize(); }

std::vector<size_t> Grid::getAlgorithmicDimensions() { return storage.getAlgorithmicDimensions(); }

void Grid::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
  this->storage.setAlgorithmicDimensions(newAlgoDims);
}

GridType Grid::stringToGridType(const std::string& gridType) {
  if (gridType.compare("linear") == 0) {
    return sgpp::base::GridType::Linear;
  } else if (gridType.compare("linearStretched") == 0) {
    return sgpp::base::GridType::LinearStretched;
  } else if (gridType.compare("linearL0Boundary") == 0) {
    return sgpp::base::GridType::LinearL0Boundary;
  } else if (gridType.compare("linearBoundary") == 0) {
    return sgpp::base::GridType::LinearBoundary;
  } else if (gridType.compare("anovaBoundary") == 0) {
    return sgpp::base::GridType::AnovaLinearBoundary;
  } else if (gridType.compare("linearStretchedBoundary") == 0) {
    return sgpp::base::GridType::LinearStretchedBoundary;
  } else if (gridType.compare("linearTruncatedBoundary") == 0) {
    return sgpp::base::GridType::LinearTruncatedBoundary;
  } else if (gridType.compare("modlinear") == 0) {
    return sgpp::base::GridType::ModLinear;
  } else if (gridType.compare("modLinearClenshawCurtis") == 0) {
    return sgpp::base::GridType::ModLinearClenshawCurtis;
  } else if (gridType.compare("poly") == 0) {
    return sgpp::base::GridType::Poly;
  } else if (gridType.compare("polyBoundary") == 0) {
    return sgpp::base::GridType::PolyBoundary;
  } else if (gridType.compare("modpoly") == 0) {
    return sgpp::base::GridType::ModPoly;
  } else if (gridType.compare("polyClenshawCurtis") == 0) {
    return sgpp::base::GridType::PolyClenshawCurtis;
  } else if (gridType.compare("polyClenshawCurtisBoundary") == 0) {
    return sgpp::base::GridType::PolyClenshawCurtisBoundary;
  } else if (gridType.compare("modPolyClenshawCurtis") == 0) {
    return sgpp::base::GridType::ModPolyClenshawCurtis;
  } else if (gridType.compare("modWavelet") == 0) {
    return sgpp::base::GridType::ModWavelet;
  } else if (gridType.compare("modBspline") == 0) {
    return sgpp::base::GridType::ModBspline;
  } else if (gridType.compare("prewavelet") == 0) {
    return sgpp::base::GridType::Prewavelet;
  } else if (gridType.compare("squareRoot") == 0) {
    return sgpp::base::GridType::SquareRoot;
  } else if (gridType.compare("periodic") == 0) {
    return sgpp::base::GridType::Periodic;
  } else if (gridType.compare("linearClenshawCurtis") == 0) {
    return sgpp::base::GridType::LinearClenshawCurtis;
  } else if (gridType.compare("linearClenshawCurtisBoundary") == 0) {
    return sgpp::base::GridType::LinearClenshawCurtisBoundary;
  } else if (gridType.compare("bspline") == 0) {
    return sgpp::base::GridType::Bspline;
  } else if (gridType.compare("bsplineBoundary") == 0) {
    return sgpp::base::GridType::BsplineBoundary;
  } else if (gridType.compare("bsplineClenshawCurtis") == 0) {
    return sgpp::base::GridType::BsplineClenshawCurtis;
  } else if (gridType.compare("wavelet") == 0) {
    return sgpp::base::GridType::Wavelet;
  } else if (gridType.compare("waveletBoundary") == 0) {
    return sgpp::base::GridType::WaveletBoundary;
  } else if (gridType.compare("fundamentalSpline") == 0) {
    return sgpp::base::GridType::FundamentalSpline;
  } else if (gridType.compare("modFundamentalSpline") == 0) {
    return sgpp::base::GridType::ModFundamentalSpline;
  } else if (gridType.compare("modBsplineClenshawCurtis") == 0) {
    return sgpp::base::GridType::ModBsplineClenshawCurtis;
  } else if (gridType.compare("linearstencil") == 0) {
    return sgpp::base::GridType::LinearStencil;
  } else if (gridType.compare("modlinearstencil") == 0) {
    return sgpp::base::GridType::ModLinearStencil;
  } else {
    std::stringstream errorString;
    errorString << "grid type '" << gridType << "' is unknown" << std::endl;
    throw base::application_exception(errorString.str().c_str());
  }
}

}  // namespace base
}  // namespace sgpp
