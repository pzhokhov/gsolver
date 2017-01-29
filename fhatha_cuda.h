#ifndef __FHATHA_CUDA_H
#define __FHATHA_CUDA_H
#include "fhatha.h"

__host__ void fhatha_runmany_cuda(fhatha_plan* p, f_complex* data, int nN, int stride, int dist);
__device __ void fhatha_cuda_kernel();

#endif