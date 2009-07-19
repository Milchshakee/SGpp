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

#include "solver/cg/ConjugateGradients.hpp"

namespace sg
{

ConjugateGradients::ConjugateGradients(size_t imax, double epsilon) : SGSolver(imax, epsilon)
{
}

ConjugateGradients::~ConjugateGradients()
{
}

void ConjugateGradients::solve(OperationMatrix& SystemMatrix, DataVector& alpha, DataVector& b, bool reuse, bool verbose, double max_threshold)
{
	this->starting();

	if (verbose == true)
	{
		std::cout << "Starting Conjugated Gradients" << std::endl;
	}

	// needed for residuum calculation
	double epsilonSquared = this->myEpsilon*this->myEpsilon;
	// number off current iterations
	this->nIterations = 0;

	// define temporal vectors
	DataVector temp(alpha.getSize());
	DataVector q(alpha.getSize());
	DataVector r(b);

	double delta_0 = 0.0;
	double delta_old = 0.0;
	double delta_new = 0.0;
	double beta = 0.0;
	double a = 0.0;

	if (verbose == true)
	{
		std::cout << "All temp variables used in CG have been intialized" << std::endl;
	}

	if (reuse == true)
	{
		q.setAll(0.0);
		SystemMatrix.mult(q, temp);
		r.sub(temp);
		delta_0 = r.dotProduct(r)*epsilonSquared;
	}
	else
	{
		alpha.setAll(0.0);
	}

	// calculate the starting residuum
	SystemMatrix.mult(alpha, temp);
	r.sub(temp);

	DataVector d(r);

	delta_old = 0.0;
	delta_new = r.dotProduct(r);

	if (reuse == false)
	{
		delta_0 = delta_new*epsilonSquared;
	}
	
	this->residuum = (delta_0/epsilonSquared);
	this->calcStarting();

	if (verbose == true)
	{
		std::cout << "Starting norm of residuum: " << (delta_0/epsilonSquared) << std::endl;
		std::cout << "Target norm:               " << (delta_0) << std::endl;
	}

	while ((this->nIterations < this->nMaxIterations) && (delta_new > delta_0) && (delta_new > max_threshold))
	{
		// q = A*d
		SystemMatrix.mult(d, q);

		// a = d_new / d.q
		a = delta_new/d.dotProduct(q);

		// x = x + a*d
		alpha.axpy(a, d);

		// Why ????
		if ((this->nIterations % 50) == 0)
		{
			// r = b - A*x
			SystemMatrix.mult(alpha, temp);
			r.copyFrom(b);
			r.sub(temp);
		}
		else
		{
			// r = r - a*q
			r.axpy(-a, q);
		}


		// calculate new deltas and determine beta
		delta_old = delta_new;
		delta_new = r.dotProduct(r);
		beta = delta_new/delta_old;
		
		this->residuum = delta_new;
		this->iterationComplete();

		if (verbose == true)
		{
			std::cout << "delta: " << delta_new << std::endl;
		}

		d.mult(beta);
		d.add(r);

		this->nIterations++;
	}

	this->residuum = delta_new;
	this->complete();

	if (verbose == true)
	{
		std::cout << "Number of iterations: " << this->nIterations << " (max. " << this->nMaxIterations << ")" << std::endl;
		std::cout << "Final norm of residuum: " << delta_new << std::endl;
	}
}

void ConjugateGradients::starting()
{
}

void ConjugateGradients::calcStarting()
{
}

void ConjugateGradients::iterationComplete()
{
}

void ConjugateGradients::complete()
{
}

}
