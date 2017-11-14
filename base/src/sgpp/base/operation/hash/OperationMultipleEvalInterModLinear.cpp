// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalInterModLinear.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

void OperationMultipleEvalInterModLinear::mult(DataVector& alpha, DataVector& result) {
  /*
  AlgorithmDGEMV<SLinearModifiedBase> op;
  LinearModifiedBasis<unsigned int, unsigned int> base;

  op.mult(storage, base, alpha, this->dataset, result);
  */
  result.setAll(0.0);


  #pragma omp parallel
  {  
    DataVector line(dataset.getNcols());
    DataVector privateResult(result.getSize());
    privateResult.setAll(0.0);

    GridStorage::grid_iterator working(storage);
    LinearModifiedBasis<unsigned int, unsigned int> basis;
    size_t dimensions = storage.getDimension();

    #pragma omp for
    for(size_t j = 0; j < dataset.getNrows(); j++){
      dataset.getRow(j,line);
      //iterate over all interactions
      for(std::vector<size_t> in:interactions){
        bool pointComputed = true;
        bool lvlComplete = false;
        size_t relLvl = 0;
        level_t* level = new level_t[in.size()];
        for (size_t i = 0; i < in.size(); i++){
          level[i] = 0;
        }
        //only view affected subspaces
        do{
          pointComputed = (!lvlComplete) & pointComputed;
          lvlComplete = false;
          for(size_t d = 1; d < dimensions; d++){
            working.push(d,1,1);
          }

          for(size_t i=0; i<in.size(); i++){
            level_t lvl = static_cast<level_t> (2 + level[i]);
            index_t idx = static_cast<index_t>(std::min(1+2*floor(line[in.at(i)] * double(1<<(lvl-1))),double(1<<lvl)-1.));
            // only hash for last element
            if(i == in.size()-1) working.set(in.at(i),lvl,idx);
            else working.push(in.at(i),lvl,idx);
          }

          // hash first dimension if there was no element in the interaction
          if(in.size() == 0) working.set(0,1,1);

          size_t seq = working.seq();

          if(!storage.isInvalidSequenceNumber(seq)){
            double value = alpha[seq];
            index_t work_index;
            level_t work_level;
            for (size_t i = 0; i < in.size(); i++){
              working.get(in.at(i),work_level,work_index);
              value *= basis.eval(work_level, work_index, line[in.at(i)]);
            }
            privateResult[j] += value;
          }

          if(in.size() == 0) break;
          //update subspace lvl
          if(level[in.size()-1] == relLvl){
            relLvl++;
            lvlComplete=true;
            for (size_t i = 0; i < in.size(); i++){
              level[i] = 0;
            }
          }
          size_t levelsum;
          //find next subspace for currentlvl
          do
          {
            levelsum = 0;
            bool carry = true;
            for(size_t i=0; i < in.size(); i++){
              if(carry) level[i]++;
              if(level[i] == relLvl+1){
                level[i] = 0;
                carry = true;
              }else carry = false;
            }
            for(size_t i=0; i < in.size(); i++){
              levelsum += level[i];
            }
          } while (levelsum != relLvl);
        }while((pointComputed||!lvlComplete));
        delete[] level;
      }
    }
    #pragma omp critical
    {
      result.add(privateResult);
    }
  }
}

void OperationMultipleEvalInterModLinear::multTranspose(DataVector& source, DataVector& result) {/*
  AlgorithmDGEMV<SLinearModifiedBase> op;
  op.mult_transposed(storage, base, source, this->dataset, result);
*/

  result.setAll(0.0);


  #pragma omp parallel
  {  
    DataVector line(dataset.getNcols());
    DataVector privateResult(result.getSize());
    privateResult.setAll(0.0);

    GridStorage::grid_iterator working(storage);
    LinearModifiedBasis<unsigned int, unsigned int> basis;
    
    size_t dimensions = storage.getDimension();

    #pragma omp for
    for(size_t j = 0; j < dataset.getNrows(); j++){
      dataset.getRow(j,line);
      //iterate over all interactions
      for(std::vector<size_t> in:interactions){
        bool pointComputed = true;
        bool lvlComplete = false;
        size_t relLvl = 0;
        level_t* level = new level_t[in.size()];
        for (size_t i = 0; i < in.size(); i++){
          level[i] = 0;
        }
        //only view affected subspaces
        do{
          pointComputed = (!lvlComplete) & pointComputed;
          lvlComplete = false;
          for(size_t d = 1; d < dimensions; d++){
            working.push(d,1,1);
          }

          for(size_t i=0; i<in.size(); i++){
            level_t lvl = static_cast<level_t> (2 + level[i]);
            index_t idx = static_cast<index_t>(std::min(1+2*floor(line[in.at(i)] * double(1<<(lvl-1))),double(1<<lvl)-1.));
            // only hash for last element
            if(i == in.size()-1) working.set(in.at(i),lvl,idx);
            else working.push(in.at(i),lvl,idx);
          }

          // hash first dimension if there was no element in the interaction
          if(in.size() == 0) working.set(0,1,1);

          size_t seq = working.seq();

          if(!storage.isInvalidSequenceNumber(seq)){
            double value = source[j];
            index_t work_index;
            level_t work_level;
            for (size_t i = 0; i < in.size(); i++){
              working.get(in.at(i),work_level,work_index);
              value *= basis.eval(work_level, work_index, line[in.at(i)]);
            }
            if(j == 0 && privateResult[seq] != 0.) std::cout << seq << "\t" << relLvl << "\t" << in.size() << std::endl;
            privateResult[seq] += value;
          }

          if(in.size() == 0) break;
          //update subspace lvl
          if(level[in.size()-1] == relLvl){
            relLvl++;
            lvlComplete=true;
            for (size_t i = 0; i < in.size(); i++){
              level[i] = 0;
            }
          }
          size_t levelsum;
          //find next subspace for currentlvl
          do
          {
            levelsum = 0;
            bool carry = true;
            for(size_t i=0; i < in.size(); i++){
              if(carry) level[i]++;
              if(level[i] == relLvl+1){
                level[i] = 0;
                carry = true;
              }else carry = false;
            }
            for(size_t i=0; i < in.size(); i++){
              levelsum += level[i];
            }
          } while (levelsum != relLvl);
        }while((pointComputed||!lvlComplete));
        delete[] level;
      }
    }
    #pragma omp critical
    {
      result.add(privateResult);
    }
  }
}


}  // namespace base
}  // namespace sgpp
