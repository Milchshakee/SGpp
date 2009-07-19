/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef SGSOLVER_HPP
#define SGSOLVER_HPP

#include "operation/OperationMatrix.hpp"
#include "solver/SGSolverInfo.hpp"
#include "solver/SGSolverResult.hpp"

namespace sg
{

/**
 * Abstract class that defines a solver used in Sparse Grids
 * Applications
 */
class SGSolver
{
private:

public:
	/// Number of Iterations needed for the solve
	size_t nIterations;
	/// Number of maximum iterations for cg
	size_t nMaxIterations;
	/// residuum
	double residuum;
	/// epsilon needed in the, e.g. final error in the iterative solver
	double myEpsilon;

	/**
	 * Std-Constructor
	 *
	 * @param nMaximumIterations number of maximum executed iterations
	 * @param epsilon the final error in the iterative solver
	 */
	SGSolver(size_t nMaximumIterations, double epsilon) : nMaxIterations(nMaximumIterations), myEpsilon(epsilon)
	{
		nIterations = 0;
		residuum = 0.0;
	}

	/**
	 * Std-Destructor
	 */
	virtual ~SGSolver() { }

	/**
	 * Pure virtual Function that defines a solve method for an iterative solver
	 *
	 * @param SystemMatrix reference to an OperationMatrix Object that implements the matrix vector multiplication
	 * @parma alpha the sparse grid's coefficients which have to be determined
	 * @param b the right hand side of the system of linear equations
	 * @param reuse identifies if the alphas, stored in alpha at calling time, should be reused
	 * @param verbose prints information during execution of the solver
	 * @param max_threshold additional abort criteria for solver
	 */
	virtual void solve(OperationMatrix& SystemMatrix, DataVector& alpha, DataVector& b, bool reuse = false, bool verbose = false, double max_threshold = -1.0) = 0;


	/**
	 * function that returns the number of needed solve steps
	 *
	 * @return the number of needed solve steps of the sovler
	 */
	size_t getNumberIterations()
	{
		return nIterations;
	}

	/**
	 * function the returns the residuum (current or final), error of the solver
	 *
	 * @return the residuum
	 */
	double getResiduum()
	{
		return residuum;
	}

	/**
	 * resets the number of maximum iterations
	 *
	 * @data iterations the new number of maximum iterations
	 */
	void setMaxIterations(size_t nIterations)
	{
		nMaxIterations = nIterations;
	}

	/**
	 * resets the epsilon, that is used in the SGSolver
	 *
	 * @data eps the new value of epsilon
	 */
	void setEpsilon(double eps)
	{
		myEpsilon = eps;
	}
};

}

#endif /* SGSOLVER_HPP */
