#ifndef _EXTENDED_MATH_H
#define _EXTENDED_MATH_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex> 

#define float_type double
using std::complex;
typedef complex<float_type> f_complex;

#ifdef _WIN32 
    inline float_type drand48() {return rand()/RAND_MAX;}
    inline void srand48(unsigned int N) {srand(N);}
	inline float_type exp2(float_type a) {return exp(log(2.0)*a);}
	inline float_type nan() {float_type c=0; return 0.0/c;}
	#define NAN (nan())
    #define isnan(a) _isnan(a) 
#endif


float_type dawson(float_type x);

float_type ellipk(float_type z);
float_type ellipe(float_type z);
void   ellipke(float_type z, float_type* k, float_type* e);
float_type oddpow(float_type x, int pw);
float_type oddroot(float_type x, int pw);
#define EM_ELLIP_TOLERANCE         1e-8


#endif
