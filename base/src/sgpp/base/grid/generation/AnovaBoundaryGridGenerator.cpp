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

#include <vector>

namespace sgpp {
namespace base {

AnovaBoundaryGridGenerator::AnovaBoundaryGridGenerator(HashGridStorage& storage)
    : storage(storage){}

AnovaBoundaryGridGenerator::~AnovaBoundaryGridGenerator() = default;

void AnovaBoundaryGridGenerator::regular(size_t level) {
  HashGenerator gen;
  gen.regularWithAnovaBoundaries(this->storage, static_cast<level_t>(level));
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

void AnovaBoundaryGridGenerator::coarsen(CoarseningFunctor& func, std::vector<size_t>* removedSeq) {
}

void AnovaBoundaryGridGenerator::coarsenNFirstOnly(CoarseningFunctor& func, size_t numFirstOnly,
  std::vector<size_t>* removedSeq, size_t minIndexConsidered) {
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
