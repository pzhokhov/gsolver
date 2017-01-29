#ifndef __FHATHA_H
#define __FHATHA_H

#include "extmath.h"

class fhatha_plan
{
	size_t N;
	float_type Nf;
	f_complex* j1; 
	
	float_type x0;
	float_type alpha;
	float_type k0;

	float_type* x;
	float_type* phi_mult;
	f_complex* phi;
	fftwt_plan f_plan; 
	
	fhatha_plan();
	
	float_type optimal_alpha(int N);
public:
	fhatha_plan(size_t N, float_type alm=1, float_type Nf=-1);
	~fhatha_plan();

    void run_one(f_complex* data);
	void run_many(f_complex* data, size_t nN, size_t stride, size_t dist);
	void run_many(f_complex* data, size_t nN) {run_many(data, nN, 1, N);}

	inline float_type x_n(size_t n) {return x[n];}
	inline float_type getNf() {return Nf;}
	void setNf(float_type newNf);

	friend void fhatha_runmany_MPI      (fhatha_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, f_complex* buf, MPI_Comm comm);

#ifdef _ENABLE_CUDA
	friend void fhatha_runmany_cuda(fhatha_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, bool do_extra_transfers = false, float_type* cuda_bufs = NULL);
	friend void fhatha_runmany_MPI_cuda (fhatha_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, f_complex* buf, MPI_Comm comm);
#endif

};



#define FHATHA_ALPHA_RELTOL 1e-8


class qdht_plan
{
   size_t N; 
   float_type V;

   float_type* m1;
   float_type* m2; 
   float_type* C; 
   float_type* x;
   f_complex* buf;
 public: 
   qdht_plan(size_t iN);
  ~qdht_plan();
   void run_one(f_complex* data);
   void run_many(f_complex* data, size_t nN, size_t stride, size_t dist);
   void run_many(f_complex* data, size_t nN) {run_many(data, nN, 1, N);}
   inline float_type x_n(size_t n) {return x[n];}
   inline float_type getNf() {return V;}
   friend void qdht_runmany_MPI (qdht_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, f_complex* buf, MPI_Comm comm); 
 
  #ifdef _ENABLE_CUDA
	friend void qdht_runmany_cuda(qdht_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, bool do_extra_transfers = false, float_type* cuda_bufs = NULL);
	friend void qdht_runmany_MPI_cuda (qdht_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, f_complex* buf, MPI_Comm comm);
 #endif  
};

#endif
