#ifndef __CUDA_EXTMATH
#define __CUDA_EXTMATH

#include <math.h>
#include "extmath.h"
#include <cuda.h>
//#include <cutil_inline.h>

//__device__ inline float_type cuda_abs2(float_type* z)                {return (*z)*(*z)+(*(z+1))*(*(z+1));}
__device__ inline float_type cuda_abs2(float_type re, float_type im) {return re*re+im*im;}

//__device__ __host__ inline int floor_(float_type x) { int c=(int)x; return (c<=x)?(c):(c-1);}
#define floor_ floor

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

__device__ inline void fft_device(float_type* data, unsigned int nn, int isign)
/*Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
be an integer power of 2 (this is not checked for!).*/

{
 unsigned int n,mmax,m,j,istep,i;
 float_type wtemp,wr,wpr,wpi,wi,theta;
 float_type tempr,tempi;
 float_type rdi, rdj, idi, idj;

 data--; 

 n=nn << 1;
 j=1;
 for (i=1;i<n;i+=2)
 {
  // This is the bit-reversal section of the routine.
  if (j > i)
  {
   SWAP(data[j],data[i]);
   //Exchange the two complex numbers.
   SWAP(data[(j+1)],data[(i+1)]);
  }
  m=n >> 1;
  while (m >= 2 && j > m)
  {
   j -= m;
   m >>= 1;
  }
  j += m;
 }
 //Here begins the Danielson-Lanczos section of the routine.
 mmax=2;
 while (n > mmax)
 {
  //Outer loop executed log2 nn times.
  istep=mmax << 1;
  theta=(float_type)isign*(6.28318530717959/(float_type)mmax);
  //Initialize the trigonometric recurrence.
  wtemp=sin_p(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin_f(theta);
  wr=1.0;
  wi=0.0;
  for (m=1;m<mmax;m+=2)
  {
   for (i=m;i<=n;i+=istep)
   {
    j=i+mmax;
	rdj = data[j]; idj = data[j+1];
	rdi = data[i]; idi = data[i+1];
    tempr=wr*rdj-wi*idj;
    tempi=wr*idj+wi*rdj;
    data[j]    =rdi-tempr;
    data[(j+1)]=idi-tempi;
    data[i]     = rdi+tempr;
    data[(i+1)] = idi+tempi;
   }
   wr=(wtemp=wr)*wpr-wi*wpi+wr;
   wi=wi*wpr+wtemp*wpi+wi;
  }
  mmax=istep;
 }
 if (isign == 1) for (i=1;i<2*nn+1;i++) data[i]/=nn;
}


__device__ inline void fft_device_strided(float_type* data, unsigned int nn, int isign, int stride)
/*Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as −1.
data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
be an integer power of 2 (this is not checked for!).*/

{
 unsigned int n,mmax,m,j,istep,i;
 float_type wtemp,wr,wpr,wpi,wi,theta;
 float_type tempr,tempi;
 float_type rdi, idi, rdj, idj; 

 data-=stride; 

 n=nn << 1;
 j=1;
 for (i=1;i<n;i+=2)
 {
  // This is the bit-reversal section of the routine.
  if (j > i)
  {
   SWAP(data[j*stride],data[i*stride]);
   //Exchange the two complex numbers.
   SWAP(data[(j+1)*stride],data[(i+1)*stride]);
  }
  m=n >> 1;
  while (m >= 2 && j > m)
  {
   j -= m;
   m >>= 1;
  }
  j += m;
 }
 //Here begins the Danielson-Lanczos section of the routine.
 mmax=2;
 while (n > mmax)
 {
  //Outer loop executed log2 nn times.
  istep=mmax << 1;
  theta=(float_type)isign*(6.28318530717959/(float_type)mmax);
  //Initialize the trigonometric recurrence.
  wtemp=sin_p(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin_p(theta);
  wr=1.0;
  wi=0.0;
  for (m=1;m<mmax;m+=2)
  {
   for (i=m;i<=n;i+=istep)
   {
    j=i+mmax;
	rdj = data[j*stride]; idj = data[(j+1)*stride];
	rdi = data[i*stride]; idi = data[(i+1)*stride];
    tempr=wr*rdj-wi*idj;
    tempi=wr*idj+wi*rdj;
    data[j*stride]    =rdi-tempr;
    data[(j+1)*stride]=idi-tempi;
    data[i*stride]     = rdi+tempr;
    data[(i+1)*stride] = idi+tempi;
   }
   wr=(wtemp=wr)*wpr-wi*wpi+wr;
   wi=wi*wpr+wtemp*wpi+wi;
  }
  mmax=istep;
 }
 if (isign == 1) for (i=1;i<2*nn+1;i++) data[i*stride]/=nn;
}


__device__ inline void getpolar(float_type reX, float_type imX, float_type* ro, float_type* phi)
{
	float_type A = sqrt_p(reX*reX + imX*imX);
	(*ro) = A; (*phi) = acos_p(reX/A)*sign(imX);
}


__device__ inline void getpolar(float_type* X, float_type* ro, float_type* phi)
{
	float_type reX = X[0], imX = X[1]; 
	getpolar(reX, imX, ro, phi);
}



__device__ inline void sqrtc(float_type reX, float_type imX, float_type* reY, float_type* imY)
{
	float_type A,phi; getpolar(reX, imX, &A, &phi); 
	float_type s,c;   sincos_f(phi/2.0, &s, &c);
	float_type absy  = sqrt_p(A);
	(*reY) = absy*c; (*imY) = absy*s;
}

__device__ inline void sqrtc(float_type* X, float_type* reY, float_type* imY)
{
	float_type reX = X[0], imX = X[1]; 
	sqrtc(reX, imX, reY, imY);
}


__device__ inline float_type device_taylor_sum(float_type x, float_type* k, int N, int startn)
{
	float_type r = 0;
	float_type xn = 1; for (int j=0; j<startn; j++) xn*=x;
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}

__device__ inline float_type device_taylor_sum1(float_type x, float_type* k, int N)
{
	float_type r = 0;
	float_type xn = x;
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}


__device__ inline float_type device_taylor_sum0(float_type x, float_type* k, int N)
{
	float_type r = 0;
	float_type xn = 1;
	for (int j=0; j<N; j++) {r+=k[j]*xn; xn*=x;}
	return r;
}


__device__ inline float_type device_sqrtHO(float_type x)
{
	float_type k[20] =  {0.5, -0.125, 0.0625, -0.0390625, 0.02734375, -0.0205078125, 0.01611328125, -0.013092041015625, 0.010910034179688, \
		-0.009273529052734, 0.008008956909180, -0.007007837295532, 0.006199240684509, -0.005535036325455, 0.004981532692909, -0.004514514002949, \
		0.004116174532101, -0.003773159987759, 0.003475278936094, -0.003214633015887};

	return device_taylor_sum1(x, k, 20);
}



__device__ inline void device_taylor_sum(float_type rex, float_type imx, float_type* rer, float_type* imr, float_type* k, int N, int startn)
{
	float_type rexn = 1, imxn = 0;
	float_type t = 0;
	float_type rer_ = 0, imr_ = 0;

	for (int j=0; j<startn; j++) 
	{
		t=rexn*rex-imxn*imx;
		imxn=rexn*imx+imxn*rex; 
		rexn=t;
	}
	for (int j=0; j<N; j++)
	{
		rer_+=k[j]*rexn;
		imr_+=k[j]*imxn;
		t=rexn*rex-imxn*imx;
		imxn=rexn*imx+imxn*rex;
		rexn=t;
	}
	(*rer)=rer_;
	(*imr)=imr_;
}

__device__ inline void device_taylor_sum1(float_type rex, float_type imx, float_type* rer, float_type* imr, float_type* k, int N)
{
	float_type rexn = rex, imxn = imx;
	float_type t = 0;
	float_type rer_ = 0, imr_ = 0;

	for (int j=0; j<N; j++)
	{
		rer_+=k[j]*rexn;
		imr_+=k[j]*imxn;
		t=rexn*rex-imxn*imx;
		imxn=rexn*imx+imxn*rex;
		rexn=t;
	}
	(*rer)=rer_;
	(*imr)=imr_;
}


__device__ inline void device_taylor_sum0(float_type rex, float_type imx, float_type* rer, float_type* imr, float_type* k, int N)
{
	float_type rexn = 1, imxn = 0;
	float_type t = 0;
	float_type rer_ = 0, imr_ = 0;

	for (int j=0; j<N; j++)
	{
		rer_+=k[j]*rexn;
		imr_+=k[j]*imxn;
		t   =rexn*rex-imxn*imx;
		imxn=rexn*imx+imxn*rex;
		rexn=t;
	}
	(*rer)=rer_;
	(*imr)=imr_;
}



__device__ inline void device_sqrtHO(float_type rex, float_type imx, float_type* rer, float_type* imr)
{
	float_type k[20] =  {0.5, -0.125, 0.0625, -0.0390625, 0.02734375, -0.0205078125, 0.01611328125, -0.013092041015625, 0.010910034179688, \
		-0.009273529052734, 0.008008956909180, -0.007007837295532, 0.006199240684509, -0.005535036325455, 0.004981532692909, -0.004514514002949, \
		0.004116174532101, -0.003773159987759, 0.003475278936094, -0.003214633015887};

	device_taylor_sum1(rex, imx, rer, imr, k, 20);
}




#endif
 
