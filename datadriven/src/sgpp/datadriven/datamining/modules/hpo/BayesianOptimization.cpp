/*
 * BayesianOptimization.cpp
 *
 *  Created on: Feb 2, 2018
 *      Author: polarbart
 */

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BayesianOptimization.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <iostream>

namespace sgpp {
namespace datadriven {

BayesianOptimization::BayesianOptimization(double firstvalue)
:kernelmatrix(1,1,1), kernelinv(1,1,1), transformedOutput(1){
	transformedOutput[0]=firstvalue;
}

double BayesianOptimization::mean(base::DataVector knew){
	return knew.dotProduct(transformedOutput);
}

double BayesianOptimization::var(base::DataVector knew, double kself){
	base::DataVector tmp(knew.size());
	kernelinv.mult(knew, tmp);
  if(knew.dotProduct(tmp)<0){
    return 1;
  }
  if(knew.dotProduct(tmp) > 1){
     // std::cout << "Error: wrong variance: "<<knew.dotProduct(tmp) <<", knew:"
     //          <<knew.toString()<<", temp: "<<tmp.toString()<<std::endl;
    return 0;
  }
	return kself - knew.dotProduct(tmp);
}

double BayesianOptimization::acquisitionEI(base::DataVector knew, double kself, double bestsofar){
    return (mean(knew) - (bestsofar-5))/var(knew, kself);
}

void BayesianOptimization::updateGP(base::DataVector knew, base::DataVector y){
	size_t idx = kernelmatrix.appendRow();
	kernelmatrix.appendCol(knew);
	kernelmatrix.setRow(idx, knew);
    kernelinv.resize(idx+1,idx+1);
  transformedOutput.resize(y.size());
    optimization::FullSLE sle(kernelmatrix);
    optimization::sle_solver::Eigen solver{};
    base::DataMatrix identity(idx+1, idx+1, 0);
    for(size_t i=0; i<idx+1;i++){
        identity.set(i,i,1);
    }
   std::cout << "vor inv solve" << std::endl;
    solver.solve(sle, identity, kernelinv);
   std::cout << "nach inv solve" << std::endl;

  // std::cout << kernelmatrix.toString() << std::endl;
  // std::cout << kernelinv.toString() << std::endl;
	kernelinv.mult(y, transformedOutput);
}

} /* namespace datadriven */
} /* namespace sgpp */
