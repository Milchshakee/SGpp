/* *****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                          *
* This file is part of the SG++ project. For conditions of distribution and    *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp           *
***************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear.hpp"
#include "base/exception/operation_exception.hpp"
#include "base/tools/AlignedMemory.hpp"

#include <omp.h>

#if defined(__SSE3__) || defined(__AVX__)
#include <immintrin.h>
#endif
#if defined(__FMA4__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

namespace sg {
  namespace parallel {

    OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear::OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear(
      sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset) :
      sg::parallel::OperationMultipleEvalVectorizedSP(storage, dataset) {
      this->level_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());
      this->index_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());
      this->mask_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());
      this->offset_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());

      storage_->getLevelIndexMaskArraysForModEval(*(this->level_), *(this->index_), *(this->mask_), *(this->offset_));

      myOCLKernels = new OCLKernels();

      _tuningMult = new sg::parallel::TwoPartitionAutoTuning(dataset->getNrows(), 128, 10);
      _tuningMultTrans = new sg::parallel::TwoPartitionAutoTuning(storage->size(), 128, 10);
    }

    OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear::~OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear() {
      delete myOCLKernels;
      delete _tuningMult;
      delete _tuningMultTrans;
    }

    void OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear::rebuildLevelAndIndex() {
      delete this->level_;
      delete this->index_;
      delete this->mask_;
      delete this->offset_;

      this->level_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());
      this->index_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());
      this->mask_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());
      this->offset_ = new sg::base::DataMatrixSP(storage_->size(), storage_->dim());

      storage_->getLevelIndexMaskArraysForModEval(*(this->level_), *(this->index_), *(this->mask_), *(this->offset_));

      myOCLKernels->resetKernels();

      _tuningMultTrans->setProblemSize(storage_->size());
      _tuningMult->resetAutoTuning();
    }

    double OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear::multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result) {
      size_t source_size = source.getSize();
      size_t dims = storage_->dim();
      size_t storageSize = storage_->size();

      double gpu_time = 0.0;
      double cpu_time = 0.0;
#ifdef _OPENMP
      int num_procs = omp_get_num_procs();
#else
      int num_procs = 1;
#endif
      double* cpu_times = new double[num_procs];

      for (int i = 0; i < num_procs; i++)
        cpu_times[i] = 0.0;

      result.setAll(0.0f);

      float* ptrSource = source.getPointer();
      float* ptrData = this->dataset_->getPointer();
      float* ptrLevel = this->level_->getPointer();
      float* ptrIndex = this->index_->getPointer();
      float* ptrMask = this->mask_->getPointer();
      float* ptrOffset = this->offset_->getPointer();
      float* ptrGlobalResult = result.getPointer();

      if (this->dataset_->getNrows() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || source_size != this->dataset_->getNrows()) {
        throw sg::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

      // split result into GPU and CPU partition
      size_t gpu_partition = storageSize - _tuningMultTrans->getPartition1Size();

      //    std::cout << gpu_partition << " " << storageSize << std::endl;

      // Do on-demand transpose
      float* ptrTransData = new float[dims * source_size];

      #pragma omp parallel for

      for (size_t n = 0; n < source_size; n++) {
        for (size_t d = 0; d < dims; d++) {
          ptrTransData[(d * source_size) + n] = ptrData[(n * dims) + d];
        }
      }

      #pragma omp parallel default(shared)
      {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
#else
        int tid = 0;
        int num_threads = 1;
#endif

        if (tid == 0) {
          if (gpu_partition > 0) {
            double loc_start = omp_get_wtime();
            myOCLKernels->multTransModMaskSPOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrMask, ptrOffset, ptrGlobalResult, source_size, storageSize, dims, gpu_partition);
            gpu_time = omp_get_wtime() - loc_start;
          }

#ifdef _OPENMP
        } else {
#endif
          size_t worksize = static_cast<size_t>((storageSize - gpu_partition) / (num_threads));
          size_t myStart = static_cast<size_t>(gpu_partition + (tid - 1) * worksize);
          size_t myEnd = static_cast<size_t>(myStart + worksize);

          if (tid == num_threads - 1)
            myEnd = storageSize;

#ifdef _OPENMP
          double start = omp_get_wtime();
#else
          myTimer_->start();
#endif
          //      #pragma omp critical
          //      {
          //        std::cout << tid << " " << myStart << " " << myEnd << " " << storageSize << std::endl;
          //      }
#if defined(__SSE3__) && !defined(__AVX__)

          for (size_t j = myStart; j < myEnd; j++) {
            __m128 res = _mm_set1_ps(0.0f);

            for (size_t i = 0; i < source_size; i += 16) {
              __m128 support_0 = _mm_load_ps(&(ptrSource[i + 0]));
              __m128 support_1 = _mm_load_ps(&(ptrSource[i + 4]));
              __m128 support_2 = _mm_load_ps(&(ptrSource[i + 8]));
              __m128 support_3 = _mm_load_ps(&(ptrSource[i + 12]));

              __m128 zero = _mm_set1_ps(0.0f);

              for (size_t d = 0; d < dims; d++) {
                __m128 eval_0 = _mm_load_ps(&(ptrTransData[(d * source_size) + i]));
                __m128 eval_1 = _mm_load_ps(&(ptrTransData[(d * source_size) + i + 4]));
                __m128 eval_2 = _mm_load_ps(&(ptrTransData[(d * source_size) + i + 8]));
                __m128 eval_3 = _mm_load_ps(&(ptrTransData[(d * source_size) + i + 12]));;

                __m128 level = _mm_load1_ps(&(ptrLevel[(j * dims) + d]));
                __m128 index = _mm_load1_ps(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                eval_0 = _mm_msub_ps(eval_0, level, index);
                eval_1 = _mm_msub_ps(eval_1, level, index);
                eval_2 = _mm_msub_ps(eval_2, level, index);
                eval_3 = _mm_msub_ps(eval_3, level, index);
#else
                eval_0 = _mm_sub_ps(_mm_mul_ps(eval_0, level), index);
                eval_1 = _mm_sub_ps(_mm_mul_ps(eval_1, level), index);
                eval_2 = _mm_sub_ps(_mm_mul_ps(eval_2, level), index);
                eval_3 = _mm_sub_ps(_mm_mul_ps(eval_3, level), index);
#endif
                __m128 mask = _mm_load1_ps(&(ptrMask[(j * dims) + d]));
                __m128 offset = _mm_load1_ps(&(ptrOffset[(j * dims) + d]));

                eval_0 = _mm_or_ps(mask, eval_0);
                eval_1 = _mm_or_ps(mask, eval_1);
                eval_2 = _mm_or_ps(mask, eval_2);
                eval_3 = _mm_or_ps(mask, eval_3);

                eval_0 = _mm_add_ps(offset, eval_0);
                eval_1 = _mm_add_ps(offset, eval_1);
                eval_2 = _mm_add_ps(offset, eval_2);
                eval_3 = _mm_add_ps(offset, eval_3);

                eval_0 = _mm_max_ps(zero, eval_0);
                eval_1 = _mm_max_ps(zero, eval_1);
                eval_2 = _mm_max_ps(zero, eval_2);
                eval_3 = _mm_max_ps(zero, eval_3);

                support_0 = _mm_mul_ps(support_0, eval_0);
                support_1 = _mm_mul_ps(support_1, eval_1);
                support_2 = _mm_mul_ps(support_2, eval_2);
                support_3 = _mm_mul_ps(support_3, eval_3);
              }

              support_0 = _mm_add_ps(support_0, support_1);
              support_2 = _mm_add_ps(support_2, support_3);
              support_0 = _mm_add_ps(support_0, support_2);

              res = _mm_add_ps(res, support_0);
            }

            res = _mm_hadd_ps(res, res);
            res = _mm_hadd_ps(res, res);

            _mm_store_ss(&(ptrGlobalResult[j]), res);
          }

#endif
#if defined(__SSE3__) && defined(__AVX__)

          for (size_t j = myStart; j < myEnd; j++) {
            __m256 res = _mm256_set1_ps(0.0f);

            for (size_t i = 0; i < source_size; i += 32) {
              __m256 support_0 = _mm256_load_ps(&(ptrSource[i + 0]));
              __m256 support_1 = _mm256_load_ps(&(ptrSource[i + 8]));
              __m256 support_2 = _mm256_load_ps(&(ptrSource[i + 16]));
              __m256 support_3 = _mm256_load_ps(&(ptrSource[i + 24]));

              __m256 zero = _mm256_set1_ps(0.0f);

              for (size_t d = 0; d < dims; d++) {
                __m256 eval_0 = _mm256_load_ps(&(ptrTransData[(d * source_size) + i]));
                __m256 eval_1 = _mm256_load_ps(&(ptrTransData[(d * source_size) + i + 8]));
                __m256 eval_2 = _mm256_load_ps(&(ptrTransData[(d * source_size) + i + 16]));
                __m256 eval_3 = _mm256_load_ps(&(ptrTransData[(d * source_size) + i + 24]));;

                __m256 level = _mm256_broadcast_ss(&(ptrLevel[(j * dims) + d]));
                __m256 index = _mm256_broadcast_ss(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                eval_0 = _mm256_msub_ps(eval_0, level, index);
                eval_1 = _mm256_msub_ps(eval_1, level, index);
                eval_2 = _mm256_msub_ps(eval_2, level, index);
                eval_3 = _mm256_msub_ps(eval_3, level, index);
#else
                eval_0 = _mm256_sub_ps(_mm256_mul_ps(eval_0, level), index);
                eval_1 = _mm256_sub_ps(_mm256_mul_ps(eval_1, level), index);
                eval_2 = _mm256_sub_ps(_mm256_mul_ps(eval_2, level), index);
                eval_3 = _mm256_sub_ps(_mm256_mul_ps(eval_3, level), index);
#endif
                __m256 mask = _mm256_broadcast_ss(&(ptrMask[(j * dims) + d]));
                __m256 offset = _mm256_broadcast_ss(&(ptrOffset[(j * dims) + d]));

                eval_0 = _mm256_or_ps(mask, eval_0);
                eval_1 = _mm256_or_ps(mask, eval_1);
                eval_2 = _mm256_or_ps(mask, eval_2);
                eval_3 = _mm256_or_ps(mask, eval_3);

                eval_0 = _mm256_add_ps(offset, eval_0);
                eval_1 = _mm256_add_ps(offset, eval_1);
                eval_2 = _mm256_add_ps(offset, eval_2);
                eval_3 = _mm256_add_ps(offset, eval_3);

                eval_0 = _mm256_max_ps(zero, eval_0);
                eval_1 = _mm256_max_ps(zero, eval_1);
                eval_2 = _mm256_max_ps(zero, eval_2);
                eval_3 = _mm256_max_ps(zero, eval_3);

                support_0 = _mm256_mul_ps(support_0, eval_0);
                support_1 = _mm256_mul_ps(support_1, eval_1);
                support_2 = _mm256_mul_ps(support_2, eval_2);
                support_3 = _mm256_mul_ps(support_3, eval_3);
              }

              support_0 = _mm256_add_ps(support_0, support_1);
              support_2 = _mm256_add_ps(support_2, support_3);
              support_0 = _mm256_add_ps(support_0, support_2);

              res = _mm256_add_ps(res, support_0);
            }

            const __m256i ldStMaskSPAVX = _mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF);

            res = _mm256_hadd_ps(res, res);
            __m256 tmp = _mm256_permute2f128_ps(res, res, 0x81);
            res = _mm256_add_ps(res, tmp);
            res = _mm256_hadd_ps(res, res);

            _mm256_maskstore_ps(&(ptrGlobalResult[j]), ldStMaskSPAVX, res);
          }

#endif
#if !defined(__SSE3__) && !defined(__AVX__)

          for (size_t j = myStart; j < myEnd; j++) {
            ptrGlobalResult[j] = 0.0f;

            for (size_t i = 0; i < source_size; i++) {
              float curSupport = ptrSource[i];

              for (size_t d = 0; d < dims; d++) {
                float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(i * dims) + d])) - (ptrIndex[(j * dims) + d]);
                unsigned int maskresult = *reinterpret_cast<unsigned int*>(&eval) | *reinterpret_cast<unsigned int*>(&(ptrMask[(j * dims) + d]));
                float masking = *reinterpret_cast<float*>( &maskresult );
                float last = masking + ptrOffset[(j * dims) + d];
                float localSupport = std::max<float>(last, 0.0f);
                curSupport *= localSupport;
              }

              ptrGlobalResult[j] += curSupport;
            }
          }

#endif
#ifdef _OPENMP
          cpu_times[tid] = omp_get_wtime() - start;
#else
          cpu_times[tid] = myTimer_->stop();
#endif
        }
      }

      for (int i = 0; i < num_procs; i++) {
        if (cpu_times[i] > cpu_time)
          cpu_time = cpu_times[i];
      }

      _tuningMultTrans->setExecutionTimes(cpu_time, gpu_time);

      double time = std::max<double>(cpu_time, gpu_time);
      //cleanup
      delete[] ptrTransData;
      delete[] cpu_times;
      return time;
    }

    double OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear::multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result) {
      size_t result_size = result.getSize();
      size_t dims = storage_->dim();
      size_t storageSize = storage_->size();

      double gpu_time = 0.0;
      double cpu_time = 0.0;
#ifdef _OPENMP
      int num_procs = omp_get_num_procs();
#else
      int num_procs = 1;
#endif
      double* cpu_times = new double[num_procs];

      for (int i = 0; i < num_procs; i++)
        cpu_times[i] = 0.0;

      result.setAll(0.0f);

      float* ptrAlpha = alpha.getPointer();
      float* ptrData = this->dataset_->getPointer();
      float* ptrResult = result.getPointer();
      float* ptrLevel = this->level_->getPointer();
      float* ptrIndex = this->index_->getPointer();
      float* ptrMask = this->mask_->getPointer();
      float* ptrOffset = this->offset_->getPointer();

      if (this->dataset_->getNrows() % 128 != 0 || result_size != this->dataset_->getNrows()) {
        throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

      // split result into GPU and CPU partition
      size_t gpu_partition = result_size - _tuningMult->getPartition1Size();

      //std::cout << gpu_partition << " " << result_size << std::endl;

      #pragma omp parallel default(shared)
      {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
#else
        int tid = 0;
        int num_threads = 1;
#endif

        if (tid == 0) {
          if (gpu_partition > 0) {
            double loc_start = omp_get_wtime();
            myOCLKernels->multModMaskSPOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrMask, ptrOffset, ptrResult, result_size, storageSize, dims, gpu_partition);
            gpu_time = omp_get_wtime() - loc_start;
          }
        } else {
          size_t worksize = static_cast<size_t>((result_size - gpu_partition) / (num_threads));
#if defined(__SSE3__) && !defined(__AVX__)
          worksize = (worksize / 16) * 16;
#endif
#if defined(__SSE3__) && defined(__AVX__)
          worksize = (worksize / 32) * 32;
#endif
          size_t myStart = static_cast<size_t>(gpu_partition + (tid - 1) * worksize);
          size_t myEnd = static_cast<size_t>(myStart + worksize);

          if (tid == num_threads - 1)
            myEnd = result_size;

#ifdef _OPENMP
          double start = omp_get_wtime();
#else
          myTimer_->start();
#endif
          //      #pragma omp critical
          //      {
          //        std::cout << tid << " " << myStart << " " << myEnd << " " << result_size << std::endl;
          //      }
#if defined(__SSE3__) && !defined(__AVX__)

          for (size_t i = myStart; i < myEnd; i += 16) {

            __m128 res_0 = _mm_load_ps(&(ptrResult[i]));
            __m128 res_1 = _mm_load_ps(&(ptrResult[i + 4]));
            __m128 res_2 = _mm_load_ps(&(ptrResult[i + 8]));
            __m128 res_3 = _mm_load_ps(&(ptrResult[i + 12]));

            // Do on-demand transpose
            float* ptrTransData = new float[dims * 16];

            for (size_t n = 0; n < 16; n++) {
              for (size_t d = 0; d < dims; d++) {
                ptrTransData[(d * 16) + n] = ptrData[((i + n) * dims) + d];
              }
            }

            for (size_t j = 0; j < storageSize; j++) {
              __m128 support_0 = _mm_load1_ps(&(ptrAlpha[j]));
              __m128 support_1 = _mm_load1_ps(&(ptrAlpha[j]));
              __m128 support_2 = _mm_load1_ps(&(ptrAlpha[j]));
              __m128 support_3 = _mm_load1_ps(&(ptrAlpha[j]));

              __m128 one = _mm_set1_ps(1.0f);
              __m128 zero = _mm_set1_ps(0.0f);
              __m128 mask = _mm_set1_ps(*fmask);

              for (size_t d = 0; d < dims; d++) {
                __m128 eval_0 = _mm_load_ps(&(ptrTransData[(d * 16)]));
                __m128 eval_1 = _mm_load_ps(&(ptrTransData[(d * 16) + 4]));
                __m128 eval_2 = _mm_load_ps(&(ptrTransData[(d * 16) + 8]));
                __m128 eval_3 = _mm_load_ps(&(ptrTransData[(d * 16) + 12]));;

                __m128 level = _mm_load1_ps(&(ptrLevel[(j * dims) + d]));
                __m128 index = _mm_load1_ps(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                eval_0 = _mm_msub_ps(eval_0, level, index);
                eval_1 = _mm_msub_ps(eval_1, level, index);
                eval_2 = _mm_msub_ps(eval_2, level, index);
                eval_3 = _mm_msub_ps(eval_3, level, index);
#else
                eval_0 = _mm_sub_ps(_mm_mul_ps(eval_0, level), index);
                eval_1 = _mm_sub_ps(_mm_mul_ps(eval_1, level), index);
                eval_2 = _mm_sub_ps(_mm_mul_ps(eval_2, level), index);
                eval_3 = _mm_sub_ps(_mm_mul_ps(eval_3, level), index);
#endif
                __m128 mask = _mm_load1_ps(&(ptrMask[(j * dims) + d]));
                __m128 offset = _mm_load1_ps(&(ptrOffset[(j * dims) + d]));

                eval_0 = _mm_or_ps(mask, eval_0);
                eval_1 = _mm_or_ps(mask, eval_1);
                eval_2 = _mm_or_ps(mask, eval_2);
                eval_3 = _mm_or_ps(mask, eval_3);

                eval_0 = _mm_add_ps(offset, eval_0);
                eval_1 = _mm_add_ps(offset, eval_1);
                eval_2 = _mm_add_ps(offset, eval_2);
                eval_3 = _mm_add_ps(offset, eval_3);

                eval_0 = _mm_max_ps(zero, eval_0);
                eval_1 = _mm_max_ps(zero, eval_1);
                eval_2 = _mm_max_ps(zero, eval_2);
                eval_3 = _mm_max_ps(zero, eval_3);

                support_0 = _mm_mul_ps(support_0, eval_0);
                support_1 = _mm_mul_ps(support_1, eval_1);
                support_2 = _mm_mul_ps(support_2, eval_2);
                support_3 = _mm_mul_ps(support_3, eval_3);
              }

              res_0 = _mm_add_ps(res_0, support_0);
              res_1 = _mm_add_ps(res_1, support_1);
              res_2 = _mm_add_ps(res_2, support_2);
              res_3 = _mm_add_ps(res_3, support_3);
            }

            delete[] ptrTransData;

            _mm_store_ps(&(ptrResult[i]), res_0);
            _mm_store_ps(&(ptrResult[i + 4]), res_1);
            _mm_store_ps(&(ptrResult[i + 8]), res_2);
            _mm_store_ps(&(ptrResult[i + 12]), res_3);
          }

#endif
#if defined(__SSE3__) && defined(__AVX__)

          for (size_t i = myStart; i < myEnd; i += 32) {

            __m256 res_0 = _mm256_load_ps(&(ptrResult[i]));
            __m256 res_1 = _mm256_load_ps(&(ptrResult[i + 8]));
            __m256 res_2 = _mm256_load_ps(&(ptrResult[i + 16]));
            __m256 res_3 = _mm256_load_ps(&(ptrResult[i + 24]));

            // Do on-demand transpose
            float* ptrTransData = new float[dims * 32];

            for (size_t n = 0; n < 32; n++) {
              for (size_t d = 0; d < dims; d++) {
                ptrTransData[(d * 32) + n] = ptrData[((i + n) * dims) + d];
              }
            }

            for (size_t j = 0; j < storageSize; j++) {
              __m256 support_0 = _mm256_broadcast_ss(&(ptrAlpha[j]));
              __m256 support_1 = _mm256_broadcast_ss(&(ptrAlpha[j]));
              __m256 support_2 = _mm256_broadcast_ss(&(ptrAlpha[j]));
              __m256 support_3 = _mm256_broadcast_ss(&(ptrAlpha[j]));

              __m256 zero = _mm256_set1_ps(0.0f);

              for (size_t d = 0; d < dims; d++) {
                __m256 eval_0 = _mm256_load_ps(&(ptrTransData[(d * 32)]));
                __m256 eval_1 = _mm256_load_ps(&(ptrTransData[(d * 32) + 8]));
                __m256 eval_2 = _mm256_load_ps(&(ptrTransData[(d * 32) + 16]));
                __m256 eval_3 = _mm256_load_ps(&(ptrTransData[(d * 32) + 24]));;

                __m256 level = _mm256_broadcast_ss(&(ptrLevel[(j * dims) + d]));
                __m256 index = _mm256_broadcast_ss(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                eval_0 = _mm256_msub_ps(eval_0, level, index);
                eval_1 = _mm256_msub_ps(eval_1, level, index);
                eval_2 = _mm256_msub_ps(eval_2, level, index);
                eval_3 = _mm256_msub_ps(eval_3, level, index);
#else
                eval_0 = _mm256_sub_ps(_mm256_mul_ps(eval_0, level), index);
                eval_1 = _mm256_sub_ps(_mm256_mul_ps(eval_1, level), index);
                eval_2 = _mm256_sub_ps(_mm256_mul_ps(eval_2, level), index);
                eval_3 = _mm256_sub_ps(_mm256_mul_ps(eval_3, level), index);
#endif
                __m256 mask = _mm256_broadcast_ss(&(ptrMask[(j * dims) + d]));
                __m256 offset = _mm256_broadcast_ss(&(ptrOffset[(j * dims) + d]));

                eval_0 = _mm256_or_ps(mask, eval_0);
                eval_1 = _mm256_or_ps(mask, eval_1);
                eval_2 = _mm256_or_ps(mask, eval_2);
                eval_3 = _mm256_or_ps(mask, eval_3);

                eval_0 = _mm256_add_ps(offset, eval_0);
                eval_1 = _mm256_add_ps(offset, eval_1);
                eval_2 = _mm256_add_ps(offset, eval_2);
                eval_3 = _mm256_add_ps(offset, eval_3);

                eval_0 = _mm256_max_ps(zero, eval_0);
                eval_1 = _mm256_max_ps(zero, eval_1);
                eval_2 = _mm256_max_ps(zero, eval_2);
                eval_3 = _mm256_max_ps(zero, eval_3);

                support_0 = _mm256_mul_ps(support_0, eval_0);
                support_1 = _mm256_mul_ps(support_1, eval_1);
                support_2 = _mm256_mul_ps(support_2, eval_2);
                support_3 = _mm256_mul_ps(support_3, eval_3);
              }

              res_0 = _mm256_add_ps(res_0, support_0);
              res_1 = _mm256_add_ps(res_1, support_1);
              res_2 = _mm256_add_ps(res_2, support_2);
              res_3 = _mm256_add_ps(res_3, support_3);
            }

            delete[] ptrTransData;

            _mm256_store_ps(&(ptrResult[i]), res_0);
            _mm256_store_ps(&(ptrResult[i + 8]), res_1);
            _mm256_store_ps(&(ptrResult[i + 16]), res_2);
            _mm256_store_ps(&(ptrResult[i + 24]), res_3);
          }

#endif
#if !defined(__SSE3__) && !defined(__AVX__)

          for (size_t i = myStart; i < myEnd; i++) {
            for (size_t j = 0; j < storageSize; j++) {
              float curSupport = ptrAlpha[j];

              for (size_t d = 0; d < dims; d++) {
                float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(i * dims) + d])) - (ptrIndex[(j * dims) + d]);
                unsigned int maskresult = *reinterpret_cast<unsigned int*>(&eval) | *reinterpret_cast<unsigned int*>(&(ptrMask[(j * dims) + d]));
                float masking = *reinterpret_cast<float*>( &maskresult );
                float last = masking + ptrOffset[(j * dims) + d];
                float localSupport = std::max<double>(last, 0.0f);
                curSupport *= localSupport;
              }

              ptrResult[i] += curSupport;
            }
          }

#endif
#ifdef _OPENMP
          cpu_times[tid] = omp_get_wtime() - start;
#else
          cpu_times[tid] = myTimer_->stop();
#endif
        }
      }

      for (int i = 0; i < num_procs; i++) {
        if (cpu_times[i] > cpu_time)
          cpu_time = cpu_times[i];
      }

      _tuningMult->setExecutionTimes(cpu_time, gpu_time);

      double time = std::max<double>(cpu_time, gpu_time);
      delete[] cpu_times;
      return time;
    }

  }

}