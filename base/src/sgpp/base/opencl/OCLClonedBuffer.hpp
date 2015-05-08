/*
 * OCLMemory.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: pfandedd
 */

#include <CL/cl.h>

#include <sgpp/globaldef.hpp>

#include "OCLManager.hpp"

namespace SGPP {
  namespace base {

    class OCLClonedBuffer {
      public:
        OCLManager& manager;
        bool initialized;
        cl_mem* bufferList;
        size_t sizeofType;
        size_t elements;

      public:

        OCLClonedBuffer(OCLManager& manager);

        bool isInitialized();

        cl_mem* getBuffer(size_t deviceNumber);

        void writeToBuffer(void* hostData, size_t* offsets = nullptr);

        void readFromBuffer(void* hostData, size_t* offsets = nullptr, size_t* ranges = nullptr);

        void initializeBuffer(void* initialValues, size_t sizeofType, size_t elements);

        void freeBuffer();

    };

  }
}


