#include <stdlib.h>
#include <fftw3.h>

#define fftwp fftw

int main()
{
 fftwp_complex* A = (fftwp_complex*)malloc(32*sizeof(fftwp_complex)); 
 double* dA = (double*)A; 

 for (int i=0; i<32; i++) dA[2*i] = (i*5 % 16) + 1.0;

 free(A); 
 return 0; 
}
