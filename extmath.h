#ifndef _EXTENDED_MATH_H
#define _EXTENDED_MATH_H

#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex> 
#include <mpi.h>
#include <stdlib.h>
//#include <math_functions.h>

#ifndef _EXTMATH_DOUBLE
#ifndef _EXTMATH_SINGLE
 //default precision is double
 #define _EXTMATH_DOUBLE
#endif
#endif

#ifndef FFTW_FLAG
#define FFTW_FLAG FFTW_ESTIMATE
#endif

 #ifndef __float128
	#define __float128 double
 #endif

 #include <fftw3.h>
 #include <fftw3-mpi.h>

using std::complex;

#ifdef _EXTMATH_SINGLE
  typedef float float_type;
  #define MPI_FLOAT_TYPE MPI_FLOAT
#else
  typedef double float_type;
  #define MPI_FLOAT_TYPE MPI_DOUBLE
#endif
  
typedef complex<float_type> f_complex;
typedef complex<double>     d_complex;
typedef complex<float>		s_complex;

    
#ifdef _WIN32 
    inline float_type frand() {return ((float_type)(rand())/((float_type)(RAND_MAX)));}
    inline void srand48(unsigned int N) {srand(N);}
	inline float_type nan() {float_type c=0.0f; return 0.0f/c;}
	#define NAN (nan())
    #define isnan(a) _isnan(a) 
#else
	inline float_type frand() {return (float_type)(drand48());}
#endif

	inline float_type frand0()  {return (frand()-(float_type)0.5);}
	inline f_complex  frand0c() {return f_complex(frand0(), frand0());}
	inline f_complex frand1c(float_type level) {return (float_type)1.0+level*frand0c();}
	inline float_type frand1(float_type level)  {return (float_type)1.0+level*frand0(); }

double dawson(double x);

double ellipk(double z);
double ellipe(double z);
void   ellipke(double z, double* k, double* e);


#ifdef _EXTMATH_SINGLE

inline float_type dawson(float_type x) {return (float_type)dawson((double)x);}
inline float_type ellipk(float_type z) {return (float_type)ellipk((double)z);}
inline float_type ellipe(float_type z) {return (float_type)ellipe((double)z);}

 #define exp_f    __expf
 #define log_f    __logf
 #define sincos_f __sincosf
 #define sin_f    __sinf
#define  cos_f    __cosf
#define  pow_f    __powf

 #define exp_p    expf
 #define log_p    logf
 #define sincos_p sincos
 #define acos_p   acosf
 #define sqrt_p   sqrt
 #define sin_p    sinf
#define  cos_p    cosf
#define  max_p    fmaxf
#define  pow_p    powf


#define __floattype2uint_rz __float2uint_rz
#define __int2floattype_rn  __int2float_rn
#define __uint2floattype_rn __uint2float_rn

inline void ellipke(float_type z, float_type* k, float_type* e) { double K,E; ellipke((double)z,&K,&E); (*k)=(float_type)K; (*e)=(float_type)E;}
#else

 #define exp_f    exp
 #define log_f    log
 #define sincos_f sincos
 #define acos_f   acos
 #define sqrt_f   sqrt
 #define max_f    max
 #define sin_f    sin
 #define cos_f    cos
#define  pow_f    pow

 #define exp_p    exp
 #define log_p    log
 #define sincos_p sincos
 #define acos_p   acos
 #define sqrt_p   sqrt
 #define max_p    max
 #define sin_p    sin
 #define cos_p    cos
 #define  pow_p   pow 

 #define __floattype2uint_rz __double2uint_rz
 #define __int2floattype_rn  __int2double_rn
 #define __uint2floattype_rn __uint2double_rn

#endif


#ifdef _EXTMATH_SINGLE
 #define fftwt_complex 				fftwf_complex
 #define fftwt_plan    				fftwf_plan

 #define fftwt_plan_dft_1d			fftwf_plan_dft_1d
 #define fftwt_plan_many_dft			fftwf_plan_many_dft
 
 #define fftwt_mpi_plan_many_dft 		fftwf_mpi_plan_many_dft
 #define fftwt_mpi_local_size_many_transposed  	fftwf_mpi_local_size_many_transposed
 
 #define fftwt_execute       			fftwf_execute
 #define fftwt_execute_dft  			fftwf_execute_dft
 #define fftwt_mpi_execute_dft			fftwf_mpi_execute_dft

 #define fftwt_destroy_plan		        fftwf_destroy_plan
 #define fftwt_malloc				fftwf_malloc

 #define fftwt_mpi_init				fftwf_mpi_init
 
#endif 

#ifdef _EXTMATH_DOUBLE
 #define fftwt_complex 				fftw_complex
 #define fftwt_plan    				fftw_plan

 #define fftwt_plan_dft_1d			fftw_plan_dft_1d
 #define fftwt_plan_many_dft			fftw_plan_many_dft
 
 #define fftwt_mpi_plan_many_dft 		fftw_mpi_plan_many_dft
 #define fftwt_mpi_local_size_many_transposed  	fftw_mpi_local_size_many_transposed
 
 #define fftwt_execute       			fftw_execute
 #define fftwt_execute_dft  			fftw_execute_dft
 #define fftwt_mpi_execute_dft			fftw_mpi_execute_dft

 #define fftwt_destroy_plan		        fftw_destroy_plan
 #define fftwt_malloc				fftw_malloc

 #define fftwt_mpi_init				fftw_mpi_init
 
#endif 



inline void* malloc_ch(size_t s) 
{
	void* p = fftwt_malloc(s); 
	if (p == NULL) 
	{
		printf("\n malloc_ch : Not enough memory. %ld bytes requested", (size_t)s); fflush(stdout);
		throw "Not enough memory. Try division into more processes or use more rough net; increasing amount of virtual memory may also help.";
	}
	return p;
}



inline float_type taylor_sum(float_type x, float_type* k, int N, int startn)
{
	float_type r = 0;
	float_type xn = 1; for (int j=0; j<startn; j++) xn*=x;
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}

inline float_type taylor_sum1(float_type x, float_type* k, int N)
{
	float_type r = 0;
	float_type xn = x;
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}


inline float_type taylor_sum0(float_type x, float_type* k, int N)
{
	float_type r = 0;
	float_type xn = 1;
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}


inline f_complex taylor_sum(f_complex x, float_type* k, int N, int startn)
{
	f_complex r = 0;
	f_complex xn = f_complex(1,0); for (int j=0; j<startn; j++) xn*=x;
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}

inline f_complex taylor_sum1(f_complex x, float_type* k, int N)
{
	f_complex r = 0;
	f_complex xn = x;
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}


inline f_complex taylor_sum0(f_complex x, float_type* k, int N)
{
	f_complex r = 0;
	f_complex xn = f_complex(1,0);
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}


inline f_complex sqrtHO(f_complex x)
{
	float_type k[20] =  {0.5, -0.125, 0.0625, -0.0390625, 0.02734375, -0.0205078125, 0.01611328125, -0.013092041015625, 0.010910034179688, \
		-0.009273529052734, 0.008008956909180, -0.007007837295532, 0.006199240684509, -0.005535036325455, 0.004981532692909, -0.004514514002949, \
		0.004116174532101, -0.003773159987759, 0.003475278936094, -0.003214633015887};

	return taylor_sum1(x, k, 20);
}



inline float_type sqrtHO(float_type x)
{
	float_type k[20] =  {0.5, -0.125, 0.0625, -0.0390625, 0.02734375, -0.0205078125, 0.01611328125, -0.013092041015625, 0.010910034179688, \
		-0.009273529052734, 0.008008956909180, -0.007007837295532, 0.006199240684509, -0.005535036325455, 0.004981532692909, -0.004514514002949, \
		0.004116174532101, -0.003773159987759, 0.003475278936094, -0.003214633015887};

	return taylor_sum1(x, k, 20);
}

	
#define sign(x)  (((x)>=0)?(1):(-1)) 

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

float_type oddpow(float_type x, int pw);
float_type oddroot(float_type x, int pw);


#define EM_ELLIP_TOLERANCE         1e-8



void least_squares(float_type* x, float_type* y, float_type* out, int N, int nN, int xstride, int xdist, int ystride, int ydist, int ostride, int odist);

inline void dump_array(void* A, size_t N, size_t size1, const char* name) {FILE* fid=fopen(name, "wb"); fwrite(A, size1, N, fid); fclose(fid);}
#define dumpar(A, N) dump_array(A, N, sizeof(A[0]), #A)


void expandtwice(f_complex* A, int stride, int dist, int Nin, int Nout, int centerpoint);
void contracttwice(f_complex* A, int stride, int dist, int Nin, int Nout, int centerpoint);



#endif
	
