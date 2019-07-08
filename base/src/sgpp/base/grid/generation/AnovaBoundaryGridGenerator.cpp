// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp>

#include <sgpp/base/grid/generation/AnovaBoundaryGridGenerator.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {

AnovaBoundaryGridGenerator::AnovaBoundaryGridGenerator(
    const AnovaHelper::AnovaComponentVector& components, HashGridStorage& storage)
    : storage(storage), components(components) {}

AnovaBoundaryGridGenerator::~AnovaBoundaryGridGenerator() = default;

void AnovaBoundaryGridGenerator::regular(size_t level) {
  HashGenerator gen;
  gen.regularWithAnovaBoundaries(this->storage, static_cast<level_t>(level), components);
}

void AnovaBoundaryGridGenerator::cliques(size_t level, size_t clique_size) {
  throw generation_exception("Method is not implemented");
}

void AnovaBoundaryGridGenerator::full(size_t level) {
  throw generation_exception("Method is not implemented");
}

void AnovaBoundaryGridGenerator::refine(RefinementFunctor& func, std::vector<size_t>* addedPoints) {
  HashRefinementBoundaries refine;
  refine.free_refine(this->storage, func, addedPoints);
}

size_t AnovaBoundaryGridGenerator::getNumberOfRefinablePoints() {
  HashRefinementBoundaries refine;
  return refine.getNumberOfRefinablePoints(this->storage);
}

void AnovaBoundaryGridGenerator::coarsen(CoarseningFunctor& func, DataVector& alpha) {
  HashCoarsening coarsen;
  coarsen.free_coarsen(this->storage, func, alpha);
}

void AnovaBoundaryGridGenerator::coarsenNFirstOnly(CoarseningFunctor& func, DataVector& alpha,
                                              size_t numFirstOnly) {
  HashCoarsening coarsen;
  coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
}

size_t AnovaBoundaryGridGenerator::getNumberOfRemovablePoints() {
  HashCoarsening coarsen;
  return coarsen.getNumberOfRemovablePoints(this->storage);
}

void AnovaBoundaryGridGenerator::refineMaxLevel(RefinementFunctor& func, size_t maxLevel) {
  HashRefinementBoundariesMaxLevel refine;
  refine.refineToMaxLevel(this->storage, func, static_cast<level_t>(maxLevel));
}

size_t AnovaBoundaryGridGenerator::getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) {
  HashRefinementBoundariesMaxLevel refine;
  return refine.getNumberOfRefinablePointsToMaxLevel(this->storage, static_cast<level_t>(maxLevel));
}

}  // namespace base
}  // namespace sgpp