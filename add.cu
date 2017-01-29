#include <stdlib.h>
#include <stdio.h>

__global__ void cudaadd(float* cA, float* cB, float* cC);

const int N = 32;
int main()
{
 int deviceN = 0; //Number of CUDA-enabled GPUs (graphics cards)
 cudaGetDeviceCount(&deviceN); 
 if (deviceN == 0) {printf("Error! No cuda-enabled devices found!"); return 1;}
 cudaSetDevice(0);  //Set 0-th device as active. Previous versions of CUDA didn't alow use of multiple GPU in a single program (except mutlti-thread programs)
// The latest CUDA release has included this featurem but I've never used it, and don't know how it works.

 float* A = (float*)malloc(N*sizeof(float));  //GeForce GPU's (like we have here) are supposed to work much faster with single-precision floating point numbers (float type) rather then 
//with double type.
 float* B = (float*)malloc(N*sizeof(float));
 float* C = (float*)malloc(N*sizeof(float));

 for (int i=0; i<N; i++) {A[i]=i%5 + i/100.0;  B[i]=2.0 + i;}  //fill in arrays

 float* cA = NULL; 
 float* cB = NULL;
 float* cC = NULL;

 cudaMalloc(&cA, N*sizeof(float)); //Allocate memory in GPU. 
 cudaMalloc(&cB, N*sizeof(float));
 cudaMalloc(&cC, N*sizeof(float));

 cudaMemcpy(cA, A, N*sizeof(float), cudaMemcpyHostToDevice); //copy arrays A and B to GPU. 
 cudaMemcpy(cB, B, N*sizeof(float), cudaMemcpyHostToDevice);
// Take care! cA and cB point to address in GPU memory. You cannot directly write there (e.g. like cA[i] = 10.1). You MUST use cudaMemcpy

 cudaadd<<<1,N>>>(cA, cB, cC);  //call GPU procedure (or "kernel"), using 1 block with N threads in block

 cudaMemcpy(C, cC, N*sizeof(float), cudaMemcpyDeviceToHost);

 for (int i=0; i<N; i++) printf("\nA[%d]=%g,  B[%d]=%g, C[%d]=%g, should be %g", i, A[i], i, B[i], i, C[i], A[i]+B[i]);

 free(A); free(B); free(C);                      //Free host arrays
 cudaFree(cA); cudaFree(cB); cudaFree(cC);       //Free GPU arrays
 return 0;
}


__global__ void cudaadd(float* cA, float* cB, float* cC)
{
 int n = threadIdx.x;
 cC[n] = cA[n]+cB[n];
}
