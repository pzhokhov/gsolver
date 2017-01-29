#include "fhatha.h"

#include <string.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

double fhatha_alpha_equation    (double z, void* params)
{
    int N = (*((int*)params));
	return (N-1)*z+log(1-exp(-z));
}

double fhatha_alpha_equation_df (double z, void* params)
{
    int N = (*((int*)params));
	return (N-1)+1.0/(exp(z)-1.0);
}
void fhatha_alpha_equation_fdf(double z, void* params, double* f, double* df) 
{
	int N = (*((int*)params));
	(*f) =  (N-1)*z+log(1.0-exp(-z));
	(*df) = (N-1)+1.0/(exp(z)-1.0);
}


fhatha_plan :: fhatha_plan() {throw "fhatha_plan: Empty plan is not allowed!";}

fhatha_plan :: fhatha_plan(size_t N_,float_type alm, float_type Nf_) : N(N_)
{
	x        = (float_type*)malloc_ch(  N*sizeof(float_type));
	j1       = (f_complex*) malloc_ch(2*N*sizeof(f_complex));  
	phi_mult = (float_type*)malloc_ch(  N*sizeof(float_type)); 
	phi      = (f_complex*) malloc_ch(2*N*sizeof(f_complex));	

	alpha = optimal_alpha((int)N); alpha*=alm;

	x0 = (1+exp(alpha))*exp(-alpha*(int)N)/2;
	k0 = (2*exp(alpha)+exp(2*alpha))/(1+exp(alpha))/(1+exp(alpha))/(1-exp(-2*alpha));

	for (int n=0; n<N;   n++) {x[n]    = x0*exp(alpha*n); phi_mult[n]=exp(alpha*((int)(1-N+n)));} phi_mult[0]*=k0;

	if (Nf_ == -1) Nf_ = 0.5/alpha;
	setNf(Nf_);

	f_plan = fftwt_plan_dft_1d(2*N, (fftwt_complex*)phi, (fftwt_complex*)phi, FFTW_FORWARD, FFTW_FLAG);
}


void fhatha_plan :: setNf(float_type newNf)
{
	Nf = newNf;
	fftwt_plan fftpl = fftwt_plan_dft_1d(2*N, (fftwt_complex*)j1, (fftwt_complex*)j1, FFTW_BACKWARD, FFTW_FLAG); 

	for (size_t n=0; n<2*N; n++) j1[n] = gsl_sf_bessel_J1(2.0*M_PI*Nf*x0*exp(alpha*((int)(n+1-N)))); 
	fftwt_execute(fftpl); for (size_t n=0; n<2*N; n++) j1[n]/=((int)2*N);
}


fhatha_plan :: ~fhatha_plan()
{
	free(x); free(phi); free(j1); free(phi_mult);
}

void fhatha_plan :: run_one(f_complex* data)
{
	for (size_t n=0; n<N-1; n++) phi[n] = phi_mult[n]*(data[n]-data[n+1]);
	phi[N-1]=data[N-1];
	for (size_t n=N; n<2*N; n++) phi[n] = 0;
        fftwt_execute(f_plan);
	for (size_t n=0; n<2*N; n++) phi[n] *= j1[n];
	fftwt_execute(f_plan);
	for (size_t n=0; n<N; n++) data[n]=phi[n]/x[n];
}

void fhatha_plan :: run_many(f_complex* data, size_t nN, size_t stride, size_t dist)
{
	for (size_t i=0; i<nN; i++)
	{
		f_complex* d = data+i*dist;
		for (size_t n=0; n<(N-1); n++) phi[n] = phi_mult[n]*(d[n*stride] - d[(n+1)*stride]);
		phi[N-1]=d[(N-1)*stride];
		for (size_t n=N; n<2*N; n++) phi[n] = 0;

		fftwt_execute(f_plan);
		for (size_t n=0; n<2*N; n++) phi[n] *= j1[n];
		fftwt_execute(f_plan);

		for (size_t n=0; n<N; n++) d[n*stride]=phi[n]/x[n];
	}
}


float_type fhatha_plan :: optimal_alpha(int N)
{
	const gsl_root_fdfsolver_type* T = gsl_root_fdfsolver_secant;
	gsl_root_fdfsolver* slv = gsl_root_fdfsolver_alloc(T);
	gsl_function_fdf func; 
	func.f   = &fhatha_alpha_equation;
	func.df  = &fhatha_alpha_equation_df;
	func.fdf = &fhatha_alpha_equation_fdf;
	func.params = &N;

	double x0=0, x1=1;
	gsl_root_fdfsolver_set(slv, &func, x1);
	
	int status = GSL_CONTINUE; 
	 
    while (status == GSL_CONTINUE)
	{
		status = gsl_root_fdfsolver_iterate(slv);
		x0 = x1;
		x1 = gsl_root_fdfsolver_root(slv);
		status = gsl_root_test_delta(x0,x1,0, FHATHA_ALPHA_RELTOL);
	}
	

	return (float_type)x1;
}


void fhatha_runmany_MPI (fhatha_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, f_complex* buf, MPI_Comm comm)
{
	int procN, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &procN);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int myN      = p->N / procN;
	int myNstart = myrank*myN;  

	int mynN     = nN/procN; 

	for (size_t i=0; i<nN;   i++) for (size_t j=0; j<myN; j++)buf[j+i*myN] = data[i*dist+j*stride];


	MPI_Alltoall(buf,  2*myN*mynN, MPI_FLOAT_TYPE, data, 2*myN*mynN, MPI_FLOAT_TYPE, comm);

	for (size_t j=0; j<p->N; j++) for (size_t i=0; i<mynN; i++) buf[j+p->N*i] = data[j%myN + (j/myN)*myN*mynN + myN*i];

	p->run_many(buf, mynN);

	for (size_t j=0; j<p->N; j++) for (size_t i=0; i<mynN; i++) data[j%myN + (j/myN)*myN*mynN + myN*i] = buf[j+p->N*i];

	MPI_Alltoall(data, 2*myN*mynN, MPI_FLOAT_TYPE, buf,  2*myN*mynN, MPI_FLOAT_TYPE, comm);

	for (size_t i=0; i<nN;   i++) for (size_t j=0; j<myN; j++) data[i*dist+j*stride] = buf[j+i*myN];
}


#ifdef _ENABLE_CUDA
void fhatha_runmany_MPI_cuda (fhatha_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, f_complex* buf, MPI_Comm comm)
{
	int procN, myrank;
	MPI_Comm_size(comm, &procN);
	MPI_Comm_rank(comm, &myrank);

	int myN      = p->N / procN;
	int myNstart = myrank*myN;  

	int mynN     = nN/procN; 

	
	for (size_t i=0; i<nN;   i++) for (size_t j=0; j<myN; j++) buf[j+i*myN] = data[i*dist+j*stride];
	MPI_Alltoall(buf,  2*myN*mynN, MPI_FLOAT_TYPE, data, 2*myN*mynN, MPI_FLOAT_TYPE, comm);
	for (size_t j=0; j<p->N; j++) for (size_t i=0; i<mynN; i++) buf[j+p->N*i] = data[j%myN + (j/myN)*myN*mynN + myN*i];

	fhatha_runmany_cuda(p, buf, mynN, 1, p->N);

	for (size_t j=0; j<p->N; j++) for (size_t i=0; i<mynN; i++) data[j%myN + (j/myN)*myN*mynN + myN*i] = buf[j+p->N*i];
	MPI_Alltoall(data, 2*myN*mynN, MPI_FLOAT_TYPE, buf,  2*myN*mynN, MPI_FLOAT_TYPE, comm);
	for (size_t i=0; i<nN;   i++) for (size_t j=0; j<myN; j++) data[i*dist+j*stride] = buf[j+i*myN];
}

#endif

qdht_plan :: qdht_plan(size_t iN)
{
  N = iN;
  x = (float_type*)malloc_ch(sizeof(float_type)*N);
  m1 = (float_type*)malloc_ch(sizeof(float_type)*N);
  m2 = (float_type*)malloc_ch(sizeof(float_type)*N);
 

  buf = (f_complex*)malloc_ch(sizeof(f_complex)*N);  
  C = (float_type*)malloc_ch(sizeof(float_type)*N*N);
  
  float_type S = gsl_sf_bessel_zero_J0(N+1);
  V = S/2.0/M_PI;

  for (int i=0; i<N; i++) {x[i]=gsl_sf_bessel_zero_J0(i+1); m1[i]=fabs(gsl_sf_bessel_J1(x[i]));}
  
  for (int i=0; i<N; i++) for (int j=0; j<N; j++)
  {
   C[i + N*j] = 2.0/S*gsl_sf_bessel_J0(x[i]*x[j]/S)/fabs(gsl_sf_bessel_J1(x[i]))/fabs(gsl_sf_bessel_J1(x[j]));
  }  
  for (int i=0; i<N; i++) x[i]/=S; 
}

qdht_plan :: ~qdht_plan()
{
 free(x); free(m1); free(m2); free(C);
}

void qdht_plan :: run_one(f_complex* data)
{
  for (int i=0; i<N; i++) data[i]/=m1[i]; 
  for (int i=0; i<N; i++) 
  {
   buf[i]=0; 
   for (int j=0; j<N; j++) buf[i] += C[i+N*j]*data[j]; 
  }
  for (int i=0; i<N; i++) buf[i]*=m1[i];  
  memcpy(data, buf, sizeof(f_complex)*N);
}


void qdht_plan :: run_many(f_complex* data, size_t nN, size_t stride, size_t dist) 
{
  for (int n=0; n<nN; n++)
  {
   for (int i=0; i<N; i++) data[n*dist + i*stride]/=m1[i]; 
   for (int i=0; i<N; i++) 
   {
    buf[i]=0; 
    for (int j=0; j<N; j++) buf[i] += C[i+N*j]*data[j*stride + dist*n]; 
   }
   for (int i=0; i<N; i++) data[n*dist + i*stride]=buf[i]*m1[i]; 
  }   
}
  


void qdht_runmany_MPI (qdht_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, f_complex* buf, MPI_Comm comm)
{
	int procN, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &procN);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int myN      = p->N / procN;
	int myNstart = myrank*myN;  

	int mynN     = nN/procN; 

	for (size_t i=0; i<nN;   i++) for (size_t j=0; j<myN; j++)buf[j+i*myN] = data[i*dist+j*stride];


	MPI_Alltoall(buf,  2*myN*mynN, MPI_FLOAT_TYPE, data, 2*myN*mynN, MPI_FLOAT_TYPE, comm);

	for (size_t j=0; j<p->N; j++) for (size_t i=0; i<mynN; i++) buf[j+p->N*i] = data[j%myN + (j/myN)*myN*mynN + myN*i];

	p->run_many(buf, mynN);

	for (size_t j=0; j<p->N; j++) for (size_t i=0; i<mynN; i++) data[j%myN + (j/myN)*myN*mynN + myN*i] = buf[j+p->N*i];

	MPI_Alltoall(data, 2*myN*mynN, MPI_FLOAT_TYPE, buf,  2*myN*mynN, MPI_FLOAT_TYPE, comm);

	for (size_t i=0; i<nN;   i++) for (size_t j=0; j<myN; j++) data[i*dist+j*stride] = buf[j+i*myN];
}

#ifdef _ENABLE_CUDA
void qdht_runmany_MPI_cuda (qdht_plan* p, f_complex* data, size_t nN, size_t stride, size_t dist, f_complex* buf, MPI_Comm comm)
{
	int procN, myrank;
	MPI_Comm_size(comm, &procN);
	MPI_Comm_rank(comm, &myrank);

	int myN      = p->N / procN;
	int myNstart = myrank*myN;  

	int mynN     = nN/procN; 

	
	for (size_t i=0; i<nN;   i++) for (size_t j=0; j<myN; j++) buf[j+i*myN] = data[i*dist+j*stride];
	MPI_Alltoall(buf,  2*myN*mynN, MPI_FLOAT_TYPE, data, 2*myN*mynN, MPI_FLOAT_TYPE, comm);
	for (size_t j=0; j<p->N; j++) for (size_t i=0; i<mynN; i++) buf[j+p->N*i] = data[j%myN + (j/myN)*myN*mynN + myN*i];

	qdht_runmany_cuda(p, buf, mynN, 1, p->N);

	for (size_t j=0; j<p->N; j++) for (size_t i=0; i<mynN; i++) data[j%myN + (j/myN)*myN*mynN + myN*i] = buf[j+p->N*i];
	MPI_Alltoall(data, 2*myN*mynN, MPI_FLOAT_TYPE, buf,  2*myN*mynN, MPI_FLOAT_TYPE, comm);
	for (size_t i=0; i<nN;   i++) for (size_t j=0; j<myN; j++) data[i*dist+j*stride] = buf[j+i*myN];
}

#endif


