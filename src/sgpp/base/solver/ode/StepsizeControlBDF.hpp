/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#ifndef STEPSIZECONTROLBDF_HPP
#define STEPSIZECONTROLBDF_HPP

#include "base/application/ScreenOutput.hpp"
#include "base/solver/ODESolver.hpp"

namespace sg
{
namespace solver
{

/**
 * This class implements a step size control using the midpoint method and BDF2
 * for solving ordinary partial equations
 *
 * For solving the system of linear equations the
 * already implemented CG-method is used
 *
 * @version $HEAD$
 */
class StepsizeControlBDF : public ODESolver
{
private:
	/// Pointer to sg::base::ScreenOutput object
	sg::base::ScreenOutput* myScreen;

	/// epsilon for the step size control
	double myEps;

public:
	/**
	 * Std-Constructer
	 *
	 * @param nTimesteps number of maximum executed iterations
	 * @param timestepSize the size of one timestep
	 * @param eps the epsilon for the step size control
	 * @param screen possible pointer to a sg::base::ScreenOutput object
	 */
	StepsizeControlBDF(size_t nTimesteps, double timestepSize, double eps, sg::base::ScreenOutput* screen = NULL);

	/**
	 * Std-Destructor
	 */
	virtual ~StepsizeControlBDF();

	virtual void solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false);
};

}
}

#endif /* STEPSIZECONTROLBDF_HPP */