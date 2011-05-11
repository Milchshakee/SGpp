/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/datadriven/DMSystemMatrixSPVectorizedIdentity.hpp"
#include "exception/operation_exception.hpp"
#include "basis/operations_factory.hpp"
using namespace sg::base;

namespace sg
{
namespace datadriven
{

DMSystemMatrixSPVectorizedIdentity::DMSystemMatrixSPVectorizedIdentity(Grid& SparseGrid, DataMatrixSP& trainData, float lambda, std::string vecMode)
{
	// handle unsupported vector extensions
	if (vecMode != "SSE" && vecMode != "AVX" && vecMode != "OCL" && vecMode != "ArBB" && vecMode != "HYBRID_SSE_OCL")
	{
		throw new operation_exception("DMSystemMatrixVectorizedIdentity : Only SSE or AVX or OCL or ArBB or HYBRID_SSE_OCL are supported vector extensions!");
	}

	resetTimers();

	// create the operations needed in ApplyMatrix
	this->vecMode = vecMode;
	this->lamb = lambda;
	this->data = new DataMatrixSP(trainData);

	if (this->vecMode == "SSE")
	{
		this->vecWidth = 24;
	}
	else if (this->vecMode == "AVX")
	{
		this->vecWidth = 48;
	}
	else if (this->vecMode == "OCL")
	{
		this->vecWidth = 128;
	}
	else if (this->vecMode == "HYBRID_SSE_OCL")
	{
		this->vecWidth = 128;
	}
	else if (this->vecMode == "ArBB")
	{
		this->vecWidth = 16;
	}
	// should not happen because this exception should have been thrown some lines upwards!
	else
	{
		throw new operation_exception("DMSystemMatrixVectorizedIdentity : Only SSE or AVX or OCL or ArBB or HYBRID_SSE_OCL are supported vector extensions!");
	}

	numTrainingInstances = data->getNrows();

	// Assure that data has a even number of instances -> padding might be needed
	size_t remainder = data->getNrows() % this->vecWidth;
	size_t loopCount = this->vecWidth - remainder;

	if (loopCount != this->vecWidth)
	{
		DataVectorSP lastRow(data->getNcols());
		for (size_t i = 0; i < loopCount; i++)
		{
			data->getRow(data->getNrows()-1, lastRow);
			data->resize(data->getNrows()+1);
			data->setRow(data->getNrows()-1, lastRow);
		}
	}

	numPatchedTrainingInstances = data->getNrows();

	if (this->vecMode != "OCL" && this->vecMode != "ArBB" && this->vecMode != "HYBRID_SSE_OCL")
	{
		data->transpose();
	}

	this->myTimer = new SGppStopwatch();

	this->B = sg::GridOperationFactory::createOperationMultipleEvalVectorizedSP(SparseGrid, this->vecMode, this->data);
}

DMSystemMatrixSPVectorizedIdentity::~DMSystemMatrixSPVectorizedIdentity()
{
	delete this->B;
	delete this->data;
	delete this->myTimer;
}

void DMSystemMatrixSPVectorizedIdentity::mult(DataVectorSP& alpha, DataVectorSP& result)
{
	DataVectorSP temp(numPatchedTrainingInstances);

    // Operation B
	this->myTimer->start();
	this->computeTimeMult += this->B->multVectorized(alpha, temp);
    this->completeTimeMult += this->myTimer->stop();

    // patch result -> set additional entries zero
    if (numTrainingInstances != temp.getSize())
    {
    	for (size_t i = 0; i < (temp.getSize()-numTrainingInstances); i++)
    	{
    		temp.set(temp.getSize()-(i+1), 0.0f);
    	}
    }

    this->myTimer->start();
    this->computeTimeMultTrans += this->B->multTransposeVectorized(temp, result);
    this->completeTimeMultTrans += this->myTimer->stop();

    result.axpy(numTrainingInstances*this->lamb, alpha);
}

void DMSystemMatrixSPVectorizedIdentity::generateb(DataVectorSP& classes, DataVectorSP& b)
{
	DataVectorSP myClasses(classes);

	// Apply padding
	if (numPatchedTrainingInstances != myClasses.getSize())
	{
		myClasses.resizeZero(numPatchedTrainingInstances);
	}

	this->myTimer->start();
	this->computeTimeMultTrans += this->B->multTransposeVectorized(myClasses, b);
	this->completeTimeMultTrans += this->myTimer->stop();
}

void DMSystemMatrixSPVectorizedIdentity::rebuildLevelAndIndex()
{
	this->B->rebuildLevelAndIndex();
}

void DMSystemMatrixSPVectorizedIdentity::resetTimers()
{
	this->completeTimeMult = 0.0;
	this->computeTimeMult = 0.0;
	this->completeTimeMultTrans = 0.0;
	this->computeTimeMultTrans = 0.0;
}

void DMSystemMatrixSPVectorizedIdentity::getTimers(double& timeMult, double& computeMult, double& timeMultTrans, double& computeMultTrans)
{
	timeMult = this->completeTimeMult;
	computeMult = this->computeTimeMult;
	timeMultTrans = this->completeTimeMultTrans;
	computeMultTrans = this->computeTimeMultTrans;
}

}
}