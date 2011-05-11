/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLogLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp"

#include "basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiDownBBLinearStretched.hpp"
#include "basis/linearstretched/noboundary/algorithm_sweep/DPhiPhiUpBBLinearStretched.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::pde;

namespace sg
{
namespace finance
{

OperationDeltaLogLinearStretched::OperationDeltaLogLinearStretched(GridStorage* storage, DataVector& coef) : UpDownOneOpDim(storage, coef)
{
}

OperationDeltaLogLinearStretched::~OperationDeltaLogLinearStretched()
{
}

void OperationDeltaLogLinearStretched::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiUpBBLinearStretched func(this->storage);
	sweep<PhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinearStretched::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	PhiPhiDownBBLinearStretched func(this->storage);
	sweep<PhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinearStretched::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiUpBBLinearStretched func(this->storage);
	sweep<DPhiPhiUpBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationDeltaLogLinearStretched::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// dphi * phi
	DPhiPhiDownBBLinearStretched func(this->storage);
	sweep<DPhiPhiDownBBLinearStretched> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}