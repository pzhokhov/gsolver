#ifdef _ENABLE_CUDA
#include "cuda_extmath.h"
#endif 

#include "ionization.h"


void calculate_NLresponse(f_complex* input, f_complex* output);
void calculate_Hlike_response(f_complex* input, f_complex* output, int N, float_type* Zeff, float_type* alpha);

#ifdef _ENABLE_CUDA
__device__ float_type* cuda_plasma_func;
__constant__ float_type cuda_plasma_factor_re;
__constant__ float_type cuda_plasma_factor_im;

__device__ float_type* cuda_raman_func;
__device__ float_type* cuda_wavenum;
__device__ float_type* cuda_omega;
__device__ float_type* cuda_hodisp;

__constant__ __device__ int cuda_N_T;
__constant__ float_type cuda_TMIN;
__constant__ float_type cuda_TMAX;
__constant__ float_type cuda_TSTEP;

__constant__ float_type cuda_OMEGA_MAX;
__constant__ float_type cuda_OMEGA_MIN;

__constant__ float_type cuda_IONIZATION_POTENTIAL;
__constant__ float_type cuda_RECOMBINATION_TAU;
__constant__ float_type cuda_IONIZATION_POTENIAL;
#ifdef MULTI_LEVEL_IONIZATION
__device__ float_type* cuda_IONIZATION_POTENTIALS;
__device__ int cuda_IONIZATION_LEVEL_N;

#endif 
__constant__ float_type cuda_AMBIENT_CARRIER_DENSITY;
__constant__ float_type cuda_PONDEROMOTIVE_COEFFICIENT;
__constant__ float_type cuda_RAMAN_FRACTION;

__constant__ float_type cuda_N4;
__constant__ float_type cuda_N2;

#ifdef THIRD_HARMONICS
__constant__ float_type cuda_th_factor; 
#endif

__constant__ float_type cuda_AVALANCHE_CROSSSECTION;
__constant__ float_type cuda_GROUP_VELOCITY;
__constant__ float_type cuda_kxstep;
__constant__ float_type cuda_kystep;

__constant__ float_type cuda_wavenum0;
__constant__ float_type cuda_omega0;


__constant__ float_type cuda_NEUTRAL_DENSITY;
#ifdef MULTIPHOTON_IONIZATION
__constant__ int cuda_K_MPI;
__constant__ float_type cuda_BETA_MPI_LN;
#else
__constant__ float_type* cuda_IONIZATION_RATE_LN;
#endif 
__device__  float_type plasma_source_function_device(float_type reA, float_type imA, float_type ro);
__device__  float_type photoionization_function_device(float_type reA, float_type imA, float_type ro);
__device__  float_type photoionization_function_device2(float_type reA, float_type imA, float_type ro1, float_type ro2);
__device__  float_type photoabsorbtion_function_device(float_type reA, float_type imA, float_type ro);
__device__  float_type avalanche_ionization_function_device(float_type reA, float_type imA, float_type ro);
__device__  float_type recombination_function_device(float_type ro);
#ifdef MULTI_LEVEL_IONIZATION
__device__ void photoionization_functionsN_device(float_type reA, float_type imA, float_type* W, int stride);
#endif

__device__  void getpolar(float_type* X, float_type* ro, float_type* phi);

__device__  void calculate_plasmadensity_small_device(float_type* field, float_type* pro, float_type* buf);
__device__  void calculate_plasmadensity_small_device_strided       (float_type* field, float_type* pro, int stride, float_type* buf);
__device__  void calculate_plasmadensity_small_device_strided_2float(float_type* field, float_type* pro, int stride, float_type* buf); 
__device__  void calculate_plasmadensity_losses_small_device(float_type* field, float_type* pro, int stride, int rostride, float_type* loss, float_type* buf);
__device__  float_type calculate_maxplasmadensity_small_device      (float_type* field, float_type* buf); 

//__device__  void calculate_single_NLresponse_kernel        (float_type* input, float_type* output, int N_T, float_type* bufs);
__device__  void calculate_single_NLresponse_kernel_strided(float_type* input, float_type* output, int N_T, float_type* bufs);
//__device__ void calculate_single_NLresponse_kernel_strided(float_type* field, float_type* out, float_type* tempb, int stride);


__global__ void calculate_Lresponse_kernel(float_type* field, float_type* buf1, float_type* kt2_, float_type ZSTEP, size_t N_T, size_t n);
__global__ void calculate_NLresponse_kernel        (f_complex* input, f_complex* output, float_type zstep, int N_T, float_type* bufs);
__global__ void calculate_NLresponse_kernel_strided(f_complex* input, f_complex* output, float_type zstep, int N_T, float_type* bufs, float_type* maA2, float_type* maxNL2);

__global__ void calculate_plasma_2float_kernel (f_complex* input, float_type* output, f_complex* buf);
__global__ void calculate_maxplasma_kernel     (f_complex* input, float_type* output, f_complex* buf); 

#ifdef YUDIN_IVANOV_CORRECTION
__device__ float_type YI_Phi_device(float_type theta, float_type g);
#endif
 
__device__ inline float_type device_hfgaussfilter(float_type omega) {float_type domega = (omega/cuda_OMEGA_MAX-1)/ABSORBTION_LAYER_WIDTH; return (domega<0)?(float_type)1.0:exp_f(-ABSORBTION_LAYER_BETA*domega*domega);} 
__device__ inline float_type device_lfgaussfilter(float_type omega) {float_type domega = (omega*(1-ABSORBTION_LAYER_WIDTH)/cuda_OMEGA_MIN-1); return (domega>0)?(float_type)1.0:exp_f(-ABSORBTION_LAYER_BETA*domega*domega);} 
//__device__ inline float_type device_hfgaussfilter(float_type omega) {float_type domega = (2.0*omega/cuda_OMEGA_MAX-1)/0.99; return exp_p(-pow_p(domega, (float_type)100));} 
//__device__ inline float_type device_lfgaussfilter(float_type omega) {return (float_type)1.0;} 

f_complex* plasma_func_;

bool cuda_do_extra_memtransfer = true; 
bool cuda_use_pinned_memory = false;  

float_type* cuda_field; 
float_type* cuda_buf1; 
float_type* cuda_buf2;

float_type* cuda_bufs;
float_type* cuda_maxA2;
float_type* cuda_maxNL2;

 size_t cuda_device_freemem=0, cuda_device_totalmem=0;
 size_t cuda_pointN = 0, cuda_piece_size = 0; 
 size_t cuda_blocksize = 64, cuda_maxbuf_len = 0;
#endif 

void calculate_NLstep(float_type* maxI, float_type* maxNL2);
void calculate_Lstep ();


void propagator()
{

#ifndef _SILENCE
#ifdef _SHOW_EVERY_STEP
	if (ISMASTER) printf("\nCalculating nonlinear responses...");fflush(stdout);
#endif
#endif
	while (true)
	{
	 float_type max_NL2 = 0;
	 float_type max_A2   = 0;


         calculate_NLstep(&max_A2, &max_NL2);

	 float_type max_NL2_0 = max_NL2, max_A2_0 = max_A2;
	 MPI_Allreduce(&max_NL2_0, &max_NL2, 1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  
	 MPI_Allreduce(&max_A2_0,  &max_A2,  1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  

         bool repeatstep = false; 	

#ifndef _ODE_EULER
	 if (ZSTEP*ZSTEP*max_NL2 >      MAX_TOLERANCE*MAX_TOLERANCE*max_A2) {repeatstep = true; STEP_N++;}
#endif
	 ZSTEP = 0.5*MAX_TOLERANCE*sqrt(max_A2/max_NL2);
#ifdef _UNIAXIAL_FINITE_DIFFERENCE
	ZSTEP = min(ZSTEP, LIN_ZSTEP);
#endif

#ifdef _SHOW_EVERY_STEP
	 if (ISMASTER) { printf("Done. (%5.2fs)Z=%e, max_A2=%e, max_NL2=%e, Zstep=%e", MPI_Wtime()-TIME_START,CURRENT_Z, max_A2, max_NL2, ZSTEP); fflush(stdout);}
#endif
		
	 if (ZSTEP < MINSTEP_RATIO*(ZNET[n_Z]-ZNET[n_Z-1]) || !(ZSTEP > 0)) 
	 {
	  printf("[%d]: max_A2=%e, max_NL2=%e, Zstep=%e", PROCESS_RANK, max_A2, max_NL2, ZSTEP);
	  throw "Minimum Z step size reached, collapse is possible!";
	 }	 
	ZSTEP = min(ZSTEP, ZNET[n_Z]-CURRENT_Z);
	
	 if (!repeatstep)
	 {

#ifndef _UNIAXIAL_FINITE_DIFFERENCE
#ifndef NO_DIFFRACTION
 #ifdef _SHOW_EVERY_STEP
	if (ISMASTER) {printf("\nRunning Fourier transform for nonlinear response array... "); fflush(stdout); }
 #endif 
#ifdef _UNIAXIAL
	#ifdef _ENABLE_CUDA
          if (!cuda_do_extra_memtransfer) ht_run((f_complex*)cuda_buf1, (f_complex*)cuda_bufs); 
	  else
        #endif
	  ht_run(BIGBUFFER1, FIELD);
		
       	 
#else
       fftwt_execute(FFT_FWPLAN_XY); 
#endif
 #ifdef _SHOW_EVERY_STEP
	if (ISMASTER) {printf("Done"); fflush(stdout); }
 #endif 

#endif
 #ifdef _SHOW_EVERY_STEP
	if (ISMASTER) {printf("\nCalculating field change during propagation step... "); fflush(stdout); }
 #endif 
 
       calculate_Lstep(); 
 #ifdef _SHOW_EVERY_STEP
	if (ISMASTER) {printf("Done."); fflush(stdout); }
 #endif 
 

#ifndef NO_DIFFRACTION
 #ifdef _SHOW_EVERY_STEP
	if (ISMASTER) {printf("\nRunning Fourier transform for field array... "); fflush(stdout); }
 #endif 
#ifdef _UNIAXIAL
	 #ifdef _ENABLE_CUDA
          if (!cuda_do_extra_memtransfer) 
          {
           ht_run((f_complex*)cuda_field, (f_complex*)cuda_bufs);
           (cudaMemcpy(FIELD, cuda_field, cuda_piece_size, cudaMemcpyDeviceToHost));
          }
	  else
        #endif
	ht_run(FIELD, BIGBUFFER1);	
#else
       fftwt_execute(FFT_BWPLAN_XY);
       for (size_t i=0; i<N_T*MY_NX*MY_NY; i++) FIELD[i] /= N_X*N_Y;       			
#endif
 #ifdef _SHOW_EVERY_STEP
	if (ISMASTER) {printf("\nDone. "); fflush(stdout); }
 #endif 
#endif 
#else
 #ifdef _SHOW_EVERY_STEP
	if (ISMASTER) {printf("\nCalculating field change during propagation step using finite differences... "); fflush(stdout); }
 #endif 
	
	calculate_Lstep_UAFD();
 #ifdef _SHOW_EVERY_STEP
	if (ISMASTER) {printf("Done."); fflush(stdout); }
 #endif 
#endif
	if (APERTURE_N > 0) if (CURRENT_Z >= APERTURE_Z[APERTURE_N - 1])
	{
	 for (int nx=0; nx<MY_NX; nx++)
	 for (int ny=0; ny<MY_NY; ny++)
	 {
 #ifdef _UNIAXIAL
 #ifndef _UNIAXIAL_FINITE_DIFFERENCE
	   float_type R = HT_PLAN->x_n(nx+MY_NXstart)*XMAX;
 #else
	   float_type R = XMAX*(nx+MY_NXstart)/N_X;
 #endif
 #else
	   float_type x = XMIN+nx*XSTEP, y = YMIN+(ny+MY_NYstart)*YSTEP;
	   float_type R = sqrt(x*x+y*y);
 #endif	
	   if (R > APERTURE_R[APERTURE_N - 1]) for (int nw=0; nw<N_T; nw++) FIELD[nw + N_T*(nx+MY_NX*ny)]=0.0; 
	   
 	 }

#ifndef _UNIAXIAL_FINITE_DIFFERENCE		
	 memcpy(BIGBUFFER2, FIELD, MY_NX*MY_NY*N_T*sizeof(f_complex)); 
#ifdef _UNIAXIAL
	 ht_run(BIGBUFFER2, BIGBUFFER1);
#else 
	 fftwt_mpi_execute_dft(FFT_FWPLAN_XY, (fftwt_complex*)BIGBUFFER2, (fftwt_complex*)BIGBUFFER2);
#endif
#endif
	 if (ISMASTER) {printf("\n Aperture #%d of radius %g is applied", APERTURE_N, APERTURE_R[APERTURE_N-1]); fflush(stdout);}
	 APERTURE_N--; 
	}

	   return; 
	 }
	}
}

void calculate_NLstep(float_type* maxI, float_type* maxNL2)
{
#ifdef _ENABLE_CUDA
 f_complex* cuda_field_piece;
 float_type* maxA2buf, *maxNL2buf; 
 #ifdef _ODE_RK4
    int bufs_size_factor = 10;
  #endif
  #ifdef _ODE_HEUN
    int bufs_size_factor = 8;
  #endif
  #ifdef _ODE_EULER
    int bufs_size_factor = 7;
  #endif
  
   
  if (STEP_N == 0) 
  {

    int devicenum = 0; 
    (cudaGetDevice(&devicenum));
    (cudaMemGetInfo(&cuda_device_freemem, &cuda_device_totalmem));


     cuda_pointN = (int)pow(2.0,floor_(log2((double)(cuda_device_freemem/N_T/(bufs_size_factor+1)/sizeof(f_complex)))));
     if (cuda_pointN > MY_NX*MY_NY) cuda_pointN=MY_NX*MY_NY;
  
     if (cuda_blocksize > cuda_pointN) cuda_blocksize=cuda_pointN;
     cuda_piece_size = N_T*cuda_pointN*sizeof(f_complex);
     cuda_maxbuf_len = sizeof(float_type)*cuda_pointN; 


     printf("\n[%d]: Launching on CUDA device %d. It has %ld bytes of memory free. Allocating for %ld points", PROCESS_RANK, devicenum, cuda_device_freemem, cuda_pointN); fflush(stdout);
  
    if (cuda_device_freemem > N_T*MY_NX*MY_NY*(bufs_size_factor+3)*sizeof(f_complex) &&false)
    {
     cuda_do_extra_memtransfer = false;

    if (ISMASTER) printf("\n Storing all data in CUDA device memory."); 
     
    (cudaMalloc((void**)&cuda_field, cuda_piece_size)); 
    (cudaMalloc((void**)&cuda_buf1, cuda_piece_size)); 
    (cudaMalloc((void**)&cuda_buf2, cuda_piece_size));
  
    (cudaMemcpy(cuda_field, FIELD,       cuda_piece_size, cudaMemcpyHostToDevice)); 
    (cudaMemcpy(cuda_buf2,  BIGBUFFER2,  cuda_piece_size, cudaMemcpyHostToDevice)); 

	
    
    (cudaMalloc((void**)&cuda_bufs,        bufs_size_factor*cuda_piece_size));
    (cudaMalloc((void**)&cuda_maxA2,       cuda_maxbuf_len)); 
    (cudaMalloc((void**)&cuda_maxNL2,      cuda_maxbuf_len)); 
   }
   if (cuda_use_pinned_memory) 
   {	
      if (ISMASTER) printf("\n Using page-locked memory. ");
      (cudaHostRegister(FIELD, cuda_piece_size, 0)); 
      (cudaHostRegister(BIGBUFFER1, cuda_piece_size, 0)); 
      (cudaHostRegister(BIGBUFFER2, cuda_piece_size, 0)); 
   }   
  }

 maxA2buf = (float_type*)malloc_ch(sizeof(float_type)*cuda_pointN); maxNL2buf = (float_type*)malloc_ch(sizeof(float_type)*cuda_pointN); 

 if (cuda_do_extra_memtransfer)
 {
   (cudaMalloc((void**)&cuda_field_piece,                  cuda_piece_size));
   (cudaMalloc((void**)&cuda_bufs,        bufs_size_factor*cuda_piece_size));
 
    (cudaMalloc((void**)&cuda_maxA2,       cuda_maxbuf_len)); 
    (cudaMalloc((void**)&cuda_maxNL2,      cuda_maxbuf_len)); 
 
 for (long i=0; i<MY_NX*MY_NY; i+= cuda_pointN)
 {
  (cudaMemcpy(cuda_field_piece, FIELD+i*N_T,      cuda_piece_size, cudaMemcpyHostToDevice)) ;
  calculate_NLresponse_kernel_strided<<< cuda_pointN/cuda_blocksize,cuda_blocksize >>>(cuda_field_piece, cuda_field_piece, ZSTEP, N_T, cuda_bufs, cuda_maxA2, cuda_maxNL2);
  (cudaMemcpy(BIGBUFFER1+i*N_T, cuda_field_piece, cuda_piece_size, cudaMemcpyDeviceToHost)) ;
  (cudaMemcpy(maxA2buf,  cuda_maxA2,  cuda_maxbuf_len, cudaMemcpyDeviceToHost)); 
  (cudaMemcpy(maxNL2buf, cuda_maxNL2, cuda_maxbuf_len, cudaMemcpyDeviceToHost)); 
  //for (int np=0; np<cuda_pointN; np++) { (*maxI) = max((*maxI), maxA2buf[np]); (*maxNL2) = max((*maxNL2), maxNL2buf[np]);}
 }
  for (long i=0; i<MY_NX*MY_NY*N_T; i++)
  {
   (*maxI)   = max((*maxI),   abs2(FIELD[i])); 
   (*maxNL2) = max((*maxNL2), abs2(BIGBUFFER1[i])); 
  }
 }
 
 else
 {
  calculate_NLresponse_kernel_strided<<< cuda_pointN/cuda_blocksize,cuda_blocksize >>>((f_complex*)cuda_field, (f_complex*)cuda_buf1, ZSTEP, N_T, cuda_bufs, cuda_maxA2, cuda_maxNL2);
  (cudaMemcpy(maxA2buf,  cuda_maxA2,  cuda_maxbuf_len, cudaMemcpyDeviceToHost)); 
  (cudaMemcpy(maxNL2buf, cuda_maxNL2, cuda_maxbuf_len, cudaMemcpyDeviceToHost)); 
  for (int np=0; np<cuda_pointN; np++) { (*maxI) = max((*maxI), maxA2buf[np]); (*maxNL2) = max((*maxNL2), maxNL2buf[np]);}
 }
 //for (int i=0; i<MY_NX*MY_NY*N_T; i++) {(*maxI) = max((*maxI), abs2(FIELD[i])); (*maxNL2)=max((*maxNL2), abs2(BIGBUFFER1[i])); }



 if (cuda_do_extra_memtransfer)
 {
  (cudaFree(cuda_maxA2)); 
  (cudaFree(cuda_maxNL2)); 
  (cudaFree(cuda_field_piece)); 
  (cudaFree(cuda_bufs)); 
 }
 free(maxA2buf); free(maxNL2buf);
#else

  float_type max_I = 0, max_NL2 = 0;
  for (int ny=0; ny<MY_NY; ny++) for (int nx=0; nx<MY_NX; nx++)
  {
	   int ofs0 = N_T*(nx+MY_NX*ny);
#ifdef NONLINEARITY_ON
 #ifdef _ODE_EULER
	  calculate_NLresponse(FIELD+ofs0, BIGBUFFER1+ofs0);
 #endif

   #ifdef _ODE_RK4
	  calculate_NLresponse(FIELD+ofs0, NL_OUTPUT1);

	  for (int nt=0; nt<N_T; nt++) NL_INPUT[nt] = FIELD[ofs0+nt] + ZSTEP*(float_type)0.5*NL_OUTPUT1[nt];
	  calculate_NLresponse(NL_INPUT, NL_OUTPUT2);

	  for (int nt=0; nt<N_T; nt++) NL_INPUT[nt] = FIELD[ofs0+nt] + ZSTEP*(float_type)0.5*NL_OUTPUT2[nt];
	  calculate_NLresponse(NL_INPUT, NL_OUTPUT3);

	  for (int nt=0; nt<N_T; nt++) NL_INPUT[nt]= FIELD[ofs0+nt]  + ZSTEP*NL_OUTPUT3[nt];
	  calculate_NLresponse(NL_INPUT, NL_OUTPUT4); 
 #endif 

 #ifdef _ODE_HEUN
	  calculate_NLresponse(FIELD+ofs0, NL_OUTPUT1);

	  for (int nt=0; nt<N_T; nt++) NL_INPUT[nt] = FIELD[ofs0+nt] + ZSTEP*NL_OUTPUT1[nt];
	  calculate_NLresponse(NL_INPUT, NL_OUTPUT2);
 #endif

#endif

	  for (int nt=0; nt<N_T; nt++)
	  {
#ifdef NONLINEARITY_ON
 #ifdef _ODE_RK4
	   BIGBUFFER1[ofs0+nt] = (NL_OUTPUT1[nt] + (float_type)2.0*NL_OUTPUT2[nt] + (float_type)2.0*NL_OUTPUT3[nt] + NL_OUTPUT4[nt])/(float_type)6.0;
 #endif
 #ifdef _ODE_HEUN
	   BIGBUFFER1[ofs0+nt] = (NL_OUTPUT1[nt] + NL_OUTPUT2[nt])/(float_type)2.0;
 #endif
#else
       BIGBUFFER1[ofs0+nt] = 0.0;
#endif
 	   max_NL2 =  max(max_NL2,  abs2(BIGBUFFER1[ofs0+nt]));
	   max_I   =  max(max_I,    abs2(FIELD     [ofs0+nt]));
	  } 
  }
  (*maxI)=max_I; (*maxNL2)=max_NL2;
#endif 
}


#ifdef _UNIAXIAL_FINITE_DIFFERENCE

void calculate_Lstep_UAFD()
{
  f_complex j = f_complex(0,1);
  MPI_Status mpistatus; 

  for (int nw=0; nw<N_T;   nw++) 
  {
    f_complex* kappa = BIGBUFFER2      +nw*(MY_NX); 
    f_complex* khi   = BIGBUFFER2+(N_T+nw)*(MY_NX);  
    f_complex khibuf, kappabuf, fieldbuf;
    //calculate alpha, beta and gamma and, finally, khi and kappa.
    if (PROCESS_RANK > 0) 
    {
//	printf("\n[%d]:Revceiving kappa and khi with tag %d from %d...", PROCESS_RANK, nw, PROCESS_RANK-1); fflush(stdout);
	MPI_Recv(&khibuf,   2, MPI_FLOAT_TYPE, PROCESS_RANK-1, nw, MPI_COMM_WORLD, &mpistatus);
	MPI_Recv(&kappabuf, 2, MPI_FLOAT_TYPE, PROCESS_RANK-1, nw, MPI_COMM_WORLD, &mpistatus);
	MPI_Recv(&fieldbuf, 2, MPI_FLOAT_TYPE, PROCESS_RANK-1, nw, MPI_COMM_WORLD, &mpistatus);
//	printf("[%d]:Success! ", PROCESS_RANK); fflush(stdout);

	float_type rm = XMAX*(MY_NXstart-1)/N_X, r = XMAX*(MY_NXstart)/N_X, rp = XMAX*(1+MY_NXstart)/N_X;
	float_type hm = r-rm, hp=rp-r; 
	
	f_complex S = -j*ZSTEP/(real(WAVENUMBER[nw]))/(float_type)4.0;
	f_complex A = S/hm*((float_type)2.0/(rp-rm)-(float_type)1.0/(r+rm));	
        f_complex C = S/hp*((float_type)2.0/(rp-rm)+(float_type)1.0/(rp+r));
	f_complex B = S*((float_type)1.0/hm/(r+rm) - (float_type)1.0/hp/(rp+r) - (float_type)2.0/hm/hp) - (float_type)1.0; 
 
	f_complex D = -A*fieldbuf - (B+(float_type)2.0)*FIELD[nw] - C*FIELD[nw+N_T] - BIGBUFFER1[nw]*ZSTEP;
   	
	khi[0]    = (D-A*khibuf)/(A*kappabuf+B); 
        kappa[0]  = -C          /(A*kappabuf+B); 
    }    
    else {khi[0]=0; kappa[0]=1;}

    for (int nx=1; nx<MY_NX; nx++)
    {
	float_type rm = XMAX*(nx+MY_NXstart-1)/N_X, r = XMAX*(nx+MY_NXstart)/N_X, rp = XMAX*(nx+1+MY_NXstart)/N_X;
	float_type hm = r-rm, hp=rp-r; 
	size_t ofs = nw+N_T*nx;
#ifndef NO_SPACE_TIME_FOCUSING
	float_type k = real(WAVENUMBER[nw]);
#else 
	float_type k = real(WAVENUMBER0); 
#endif
	f_complex S = -j*ZSTEP/(real(WAVENUMBER[nw]))/(float_type)4.0;
	f_complex A = S/hm*((float_type)2.0/(rp-rm)-(float_type)1.0/(r+rm));	
        f_complex C = S/hp*((float_type)2.0/(rp-rm)+(float_type)1.0/(rp+r));
	f_complex B = S*((float_type)1.0/hm/(r+rm) - (float_type)1.0/hp/(rp+r) - (float_type)2.0/hm/hp) - (float_type)1.0; 
 
	f_complex D = -A*FIELD[ofs-N_T] - (B+(float_type)2.0)*FIELD[ofs] - C*FIELD[ofs+N_T] - BIGBUFFER1[ofs]*ZSTEP;
   	
	khi[nx]    = (D-A*khi[nx-1])/(A*kappa[nx-1]+B); 
        kappa[nx]  =-C              /(A*kappa[nx-1]+B);
    }
    if (PROCESS_RANK < PROCESS_N-1) 
    {
//	printf("\n[%d]:Sending kappa and khi with tag %d to %d...", PROCESS_RANK, nw, PROCESS_RANK+1); fflush(stdout);
    	MPI_Send(khi+MY_NX-1,            2, MPI_FLOAT_TYPE, PROCESS_RANK+1, nw, MPI_COMM_WORLD);
    	MPI_Send(kappa+MY_NX-1,          2, MPI_FLOAT_TYPE, PROCESS_RANK+1, nw, MPI_COMM_WORLD);
	MPI_Send(FIELD+nw+N_T*(MY_NX-1), 2, MPI_FLOAT_TYPE, PROCESS_RANK+1, nw, MPI_COMM_WORLD); 
//	printf("[%d]:Success!", PROCESS_RANK); fflush(stdout);
    }	
  }

  for (int nw=0; nw<N_T; nw++)
  {
    f_complex* kappa = BIGBUFFER2      +nw*(MY_NX);  
    f_complex* khi   = BIGBUFFER2+(N_T+nw)*(MY_NX);  
    if (PROCESS_RANK < PROCESS_N-1) MPI_Recv(FIELD+nw+MY_NX*N_T, 2, MPI_FLOAT_TYPE, PROCESS_RANK+1, nw, MPI_COMM_WORLD, &mpistatus);
    else FIELD[nw+MY_NX*N_T] = 0;//-(khi[MY_NX-1])/(kappa[MY_NX-1]-(float_type)1.0);  //Dierichlet boundary outer boundary condition

    for (int nx=MY_NX-1; nx>=0; nx--)
    {
      size_t ofs = nw+N_T*nx; 
      FIELD[ofs] = FIELD[ofs+N_T]*kappa[nx]+khi[nx];	
    }
    if (PROCESS_RANK > 0) MPI_Send(FIELD+nw, 2, MPI_FLOAT_TYPE, PROCESS_RANK-1, nw, MPI_COMM_WORLD);
  }
  
  for (int nw=0; nw<N_T; nw++) for (int nx=0; nx<=MY_NX; nx++) FIELD[nw+N_T*nx] *= exp(ZSTEP*(-j*HO_DISPERSION[nw]+imag(WAVENUMBER[nw])));
}	

#endif


void calculate_Lstep()
{
#ifdef _ENABLE_CUDA

 float_type* kt2 = (float_type*)malloc_ch(sizeof(float_type)*N_Y*MY_NX_FT);
 float_type* cuda_kt2;
 if (cuda_do_extra_memtransfer) (cudaMalloc((void**)&cuda_kt2, sizeof(float_type)*N_Y*MY_NX_FT));
 else cuda_kt2 = cuda_bufs;

#ifndef _UNIAXIAL 
 float_type kxstep = (float_type)2.0*M_PI/(XMAX-XMIN); 
 float_type kystep = (float_type)2.0*M_PI/(YMAX-YMIN); 
#endif 

 for (size_t nx=0; nx<MY_NX_FT; nx++) for (size_t ny=0; ny<N_Y; ny++)
 {
#ifdef _UNIAXIAL
   float_type kt = 2*M_PI*HT_PLAN->x_n(nx+MY_NXstart_FT)*HT_PLAN->getNf()/XMAX;
   kt2[nx+MY_NX_FT*ny] = kt*kt;
#else
  size_t nx_g = nx+MY_NXstart_FT;   
  float_type kx=0; if (nx_g<=N_X/2) kx=kxstep*nx_g; else kx=-kxstep*(N_X-nx_g);
  float_type ky=0; if (ny<=N_Y/2)   ky=kystep*ny;   else ky=-kystep*(N_Y-ny);
  kt2[ny+N_Y*nx] = kx*kx+ky*ky;
#endif 
 }
 (cudaMemcpy(cuda_kt2, kt2, sizeof(float_type)*N_Y*MY_NX_FT, cudaMemcpyHostToDevice));
 free(kt2);

 float_type* cuda_field_piece, *cuda_buf_piece; 

 if (cuda_do_extra_memtransfer)
 {
  (cudaMalloc((void**)&cuda_field_piece,    cuda_piece_size));
  (cudaMalloc((void**)&cuda_buf_piece,      cuda_piece_size));
  for (size_t i=0; i<N_Y*MY_NX_FT; i+= cuda_pointN)
  {
   (cudaMemcpy(cuda_field_piece, BIGBUFFER2+i*N_T,      cuda_piece_size, cudaMemcpyHostToDevice)) ;
   (cudaMemcpy(cuda_buf_piece,   BIGBUFFER1+i*N_T, cuda_piece_size, cudaMemcpyHostToDevice)) ;
   calculate_Lresponse_kernel<<< cuda_pointN/cuda_blocksize,cuda_blocksize >>>(cuda_field_piece, cuda_buf_piece, cuda_kt2, ZSTEP, N_T, i);
   (cudaMemcpy(FIELD+i*N_T,      cuda_field_piece, cuda_piece_size, cudaMemcpyDeviceToHost));
   (cudaMemcpy(BIGBUFFER2+i*N_T, cuda_field_piece, cuda_piece_size, cudaMemcpyDeviceToHost));
  }
 
 (cudaFree(cuda_field_piece)); 
 (cudaFree(cuda_buf_piece)); 
 (cudaFree(cuda_kt2));
 }
 else 
 {
   calculate_Lresponse_kernel<<< cuda_pointN/cuda_blocksize,cuda_blocksize >>>(cuda_buf2, cuda_buf1, cuda_kt2, ZSTEP, N_T, 0);
   (cudaMemcpy(cuda_field, cuda_buf2, cuda_piece_size, cudaMemcpyDeviceToDevice));   
 //  (cudaMemcpyAsync(BUGBUFFER2, cuda_buf2, cuda_piece_size, cudaMemcpyDeviceToHost)); 
 }
#else

#ifndef _UNIAXIAL 
 float_type kxstep = (float_type)2.0*M_PI/(XMAX-XMIN); 
 float_type kystep = (float_type)2.0*M_PI/(YMAX-YMIN); 
#endif 

 f_complex j=f_complex(0,1);
 for (size_t nx=0; nx<MY_NX_FT; nx++) for (size_t ny=0; ny<N_Y; ny++)
 {
   float_type kt2 = 0;
#ifdef _UNIAXIAL
   float_type kt = 2*M_PI*HT_PLAN->x_n(nx+MY_NXstart_FT)*HT_PLAN->getNf()/XMAX;
   kt2 = kt*kt; 
#else
  size_t nx_g = nx+MY_NXstart_FT;   
  float_type kx=0; if (nx_g<=N_X/2) kx=kxstep*nx_g; else kx=-kxstep*(N_X-nx_g);
  float_type ky=0; if (ny<=N_Y/2)   ky=kystep*ny;   else ky=-kystep*(N_Y-ny);
  kt2 = kx*kx+ky*ky;
#endif 

  for (size_t nw = 0; nw<N_T;     nw++)
  {
    size_t ofs = (nw + N_T*(ny+N_Y*nx));
	float_type w = OMEGA[nw];
#ifndef NO_SPACE_TIME_FOCUSING
	f_complex  k = WAVENUMBER[nw]; 
#else 
	f_complex  k = WAVENUMBER0;
#endif

	if (w<OMEGA_MIN || w > OMEGA_MAX)     {FIELD[ofs]=0.0;continue;}
#ifndef NO_DIFFRACTION
	if (kt2>MAX_KT2*real(k*k)) {FIELD[ofs]=0.0;continue;} 
#endif 
    //float_type f = lfgaussfilter(w)*hfgaussfilter(w);

#ifndef NO_DIFFRACTION
	float_type k1 = real(k), k2=imag(k);
#ifdef PARABOLICAL_DIFRACTION
	BIGBUFFER2[ofs] += ZSTEP*BIGBUFFER1[ofs];
	FIELD[ofs] = (BIGBUFFER2[ofs]*exp(ZSTEP*(j*(kt2/2/k1-HO_DISPERSION[nw])-k2)));// + ZSTEP*BIGBUFFER1[ofs]);
#else
    //f_complex kz = sqrt(k*k - kt2);   	     
	float_type sqHO = sqrtHO(-kt2/k1/k1);
	BIGBUFFER2[ofs] += ZSTEP*BIGBUFFER1[ofs]/(1+sqHO);
	FIELD[ofs] = (BIGBUFFER2[ofs]*exp(ZSTEP*(j*(-HO_DISPERSION[nw]-sqHO*k1)+k2)));// + ZSTEP*BIGBUFFER1[ofs]/(1+sqHO));
#endif
#else
	FIELD[ofs] = BIGBUFFER2[ofs]*exp(-ZSTEP*j*HO_DISPERSION[nw]) + ZSTEP*BIGBUFFER1[ofs];
#endif
	BIGBUFFER2[ofs]=FIELD[ofs];
  }
 }

#endif
}


void calculate_NLresponse(f_complex* input, f_complex* output)  
{
	//This function calculates nonlinear response at one spatial net point.
	f_complex j = f_complex(0,1); 

#ifndef NONLINEARITY_ON
        for (int nt=0; nt<N_T; nt++) output[nt] = 0;
        return;
#endif

#ifndef H_LIKE_RESPONSE
	fftwt_execute_dft(FFT_BWPLAN_T, (fftwt_complex*)input, (fftwt_complex*)FIELD_REAL_SMALL);
	fftwt_Nnormalize(1, FIELD_REAL_SMALL);
#ifndef NO_PLASMARESPONSE
	calculate_plasmadensity_losses_small(FIELD_REAL_SMALL, (float_type*)NL_SMALLBUFFER3, 2, NL_SMALLBUFFER1, NL_SMALLBUFFER4);           //calculate plasma density at this point
	

	//float_type tstep = (TMAX - TMIN)/N_T;

	for (int nt=0; nt<N_T; nt++)
	{
	 float_type ro = real(NL_SMALLBUFFER3[nt]);
#ifndef PLASMA_FULL_DISPERSION
	 f_complex E = FIELD_REAL_SMALL[nt];
#ifdef PLASMA_DISPERSION
	 NL_SMALLBUFFER3[nt] = ro*E;													 //multiply plasma density by field
 #else
	 f_complex Fro = ro*PLASMA_FACTOR, sFro = 0;
	 if (fabs(real(Fro))<0.1) sFro = sqrtHO(-Fro);
	 else sFro = sqrt((float_type)1.0-Fro)-(float_type)1.0;

	 NL_SMALLBUFFER1[nt] +=  - j*WAVENUMBER0*sFro*E;
 #endif 
#else
	 f_complex w0  = exp((float_type)2.0*(float_type)M_PI*j*((float_type)nt)/((float_type)N_T));
	 f_complex M   = 1.0/N_T;
	 for (int nw=0; nw<N_T; nw++)
	 { 
		f_complex Fro_ = PLASMA_FUNC[nw]*ro, sFro_ = 0;   
		if (fabs(real(Fro_))<0.1) sFro_ = sqrtHO(+Fro_);
		else sFro_ = sqrt((float_type)1.0+Fro_)-(float_type)1.0;

		NL_SMALLBUFFER1[nt]+=-j*WAVENUMBER[nw]*sFro_*input[nw]*M;
		M*=w0;
	 } 
	
#endif 
	}
	
 #ifdef PLASMA_DISPERSION
	fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)NL_SMALLBUFFER2);    //These two lines execute forward Fourier transform 1->2 and 3->1 
	fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER3, (fftwt_complex*)NL_SMALLBUFFER1);    //i.e. for ro*E and for gamma*E, where ro is plasma density, E - field and gamma - PA losses

	for (int nw=0; nw<N_T;   nw++)
	{		
	 NL_SMALLBUFFER1[nw] *= PLASMA_FUNC[nw];
	 NL_SMALLBUFFER2[nw] += NL_SMALLBUFFER1[nw];            //now NL_SMALLBUFFER2 contains all ionization responses.
	}
 #else
	fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)NL_SMALLBUFFER2);
 #endif 
#else 
	for (int nw=0; nw<N_T; nw++) {NL_SMALLBUFFER2[nw]=0; NL_SMALLBUFFER4[nw]=1.0;}
#endif
	
    if (NONLIN_REFRINDEX != 0)
    {
	 for (int nt=0; nt<N_T;  nt++) NL_SMALLBUFFER3[nt]  = abs2(FIELD_REAL_SMALL[nt]);  //put intensity into NL_SMALLBUFFER2 
	 if (RAMAN_FRACTION > 0.001)					     
	 {
	  //Calculate delayed nonlinearity response using intensity Fourier-transform
	  fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER3, (fftwt_complex*)NL_SMALLBUFFER1);
	  for (int nw=0; nw<N_T; nw++) NL_SMALLBUFFER1[nw] *= RAMAN_FUNCTION[nw];
	  fftwt_execute_dft(FFT_BWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)NL_SMALLBUFFER3);
	  fftwt_Nnormalize(1,NL_SMALLBUFFER3);
	  for (int nt=0; nt<N_T; nt++) { f_complex  E=FIELD_REAL_SMALL[nt]; NL_SMALLBUFFER3[nt]*=E;} 
	}
	else
	for (int nt=0; nt<N_T; nt++) { f_complex  E=FIELD_REAL_SMALL[nt]; NL_SMALLBUFFER3[nt]*=E;}

	fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER3, (fftwt_complex*)NL_SMALLBUFFER1); 
	for (int nw=0; nw<N_T; nw++) { NL_SMALLBUFFER2[nw] += -j*KERR_PROFILE[nw]*NL_SMALLBUFFER1[nw];}
	
#ifdef THIRD_HARMONICS
	for (int nt=0; nt<N_T; nt++) 
 	{
      	 float_type carrier_phase = (OMEGA0 == OMEGA[0])?(-2.0*OMEGA0*(TMIN+nt*(TMIN-TMAX)/N_T)):0;
 	 f_complex  E=FIELD_REAL_SMALL[nt]; NL_SMALLBUFFER3[nt]=E*E*E*exp(j*carrier_phase);
	}
	fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER3, (fftwt_complex*)NL_SMALLBUFFER1); 
	for (int nw=0; nw<N_T; nw++) { NL_SMALLBUFFER2[nw] += -j*KERR_TH_PROFILE[nw]*NL_SMALLBUFFER1[nw];}
#endif
   }
   
   if (NONLIN_REFRINDEX4 != 0)
   {
	for (int nt=0; nt<N_T; nt++) 
	{ f_complex  E=FIELD_REAL_SMALL[nt]; float_type I = abs2(E);
          NL_SMALLBUFFER3[nt]=I*I*E;
#ifdef THIRD_HARMONICS 
      	 float_type carrier_phase = (OMEGA0 == OMEGA[0])?(-2.0*OMEGA0*(TMIN+nt*(TMIN-TMAX)/N_T)):0;
          NL_SMALLBUFFER3[nt] += TH_FACTOR*(I*E*E*E*exp(j*carrier_phase)/2.0 + E*E*E*E*E*exp(2.0*j*carrier_phase)); 
#endif
        }
	fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER3, (fftwt_complex*)NL_SMALLBUFFER1); 
#ifndef NO_SHOCK
	for (int nw=0; nw<N_T; nw++) { NL_SMALLBUFFER2[nw] += -j*OMEGA[nw]/LIGHT_VELOCITY*NL_SMALLBUFFER3[nw];}		
#else
	for (int nw=0; nw<N_T; nw++) { NL_SMALLBUFFER2[nw] += -j*OMEGA[0]/LIGHT_VELOCITY*NL_SMALLBUFFER3[nw];}		
	
#endif
   }

/*
	 for (int nt=0; nt<N_T;  nt++) 
	 {
	  f_complex E = FIELD_REAL_SMALL[nt];
	  float_type I = abs2(E);
#ifndef THIRD_HARMONICS
	  NL_SMALLBUFFER3[nt] *= E*NL_SMALLBUFFER4[nt];
	  NL_SMALLBUFFER3[nt] += NONLIN_REFRINDEX4*I*I*E*NL_SMALLBUFFER3[nt];

#else
      	  float_type carrier_phase = (OMEGA0 == OMEGA[0])?(-2.0*OMEGA0*(TMIN+nt*(TMIN-TMAX)/N_T)):0;
	  if (NL_SMALLBUFFER3[nt]) = NL

	  NL_SMALLBUFFER3[nt] = NONLIN_REFRINDEX*(NL_SMALLBUFFER3[nt] + (1-RAMAN_FRACTION)*((float_type)(1.0/3.0))*TH_FACTOR*E*E*exp(j*carrier_phase))*E*NL_SMALLBUFFER4[nt];
	  NL_SMALLBUFFER3[nt] += NONLIN_REFRINDEX4*(I*I + TH_FACTOR*((float_type)(1.0/2.0)*exp(j*carrier_phase)*I*E*E + (float_type)(1.0/10.0)*exp((float_type)2.0*j*carrier_phase)*E*E*E*E))*E;
#endif
	 }

	 fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER3, (fftwt_complex*)NL_SMALLBUFFER1);

	 for (int nw=0; nw<N_T; nw++) 
	 {
#ifndef NO_SHOCK
	    float_type w = OMEGA[nw];
#else
	    float_type w = OMEGA0;
#endif
	    output[nw]	= -j*(KERR_PROFILE*NL_SMALLBUFFER1[nw]) + NL_SMALLBUFFER2[nw];
         }

    }
    else*/
   for (int nw=0; nw<N_T; nw++) output[nw] = NL_SMALLBUFFER2[nw];
#else
   float_type Zeff[2] = {1, 1}; 
   float_type alpha[2] = {1, 1}; 

   calculate_Hlike_response(input, output, 2, Zeff, alpha); 
#endif
 //  for (int nw=0; nw<N_T; nw++) {float_type w_ = OMEGA[nw]; output[nw]*=hfgaussfilter(w_)*lfgaussfilter(w_);}

}


void calculate_Hlike_response(f_complex* input, f_complex* output, int N, float_type* Zeff, float_type* alpha)
{
 float_type Xmax = 20; 
 int Nx   = 256; 
 float_type dt0 = 0.2; 
 f_complex j = f_complex(0,1); 

 int Nt   = (TMAX-TMIN)/ATOMIC_TIME/dt0; 
 int tfac = (int)exp2(ceil(log2((float_type)Nt/(float_type)N_T))); 
 Nt = N_T*tfac; 

 float_type dt =(TMAX-TMIN)/Nt/ATOMIC_TIME; 
 
 f_complex* E  = (f_complex*)malloc_ch(Nt*sizeof(f_complex)); 
 f_complex* Er = (f_complex*)malloc_ch(Nt*sizeof(f_complex)); 
 
 memcpy(E,            input,        N_T/2*sizeof(f_complex)); 
 memcpy(E+(Nt-N_T/2), input+N_T/2,  N_T/2*sizeof(f_complex)); 
 for (int nt=N_T/2; nt<(Nt-N_T/2); nt++) E[nt]=0; 

 fftwt_plan plan = fftwt_plan_dft_1d(Nt, (fftwt_complex*)E, (fftwt_complex*)Er, FFTW_BACKWARD, FFTW_ESTIMATE); 
 fftwt_execute(plan); 
 fftwt_destroy_plan(plan); 

 free(E); 

#ifndef _UNWRAP_FRQUENCIES
 for (int nt=0; nt<Nt; nt++) Er[nt] *= exp(j*OMEGA0*(TMIN+dt*nt*ATOMIC_TIME))*(float_type)FIELD_DENOM/(ATOMIC_FIELD_*N_T);
#else
 for (int nt=0; nt<Nt; nt++) Er[nt] *= (float_type)FIELD_DENOM/(ATOMIC_FIELD_);
#endif

 float_type dx = Xmax/Nx; 

 f_complex* f = (f_complex*)malloc_ch(2*N*Nx*sizeof(f_complex));   for (int n=0; n<2*N; n++) for (int nx=0; nx<Nx; nx++) f[nx+n*Nx]=1.0; 
 
 f_complex* gamma = (f_complex*)malloc_ch(2*N*Nx*sizeof(f_complex));
 f_complex* kappa = (f_complex*)malloc_ch(2*N*Nx*sizeof(f_complex)); 

 float_type* Ixn  = (float_type*)malloc_ch(6*N*sizeof(float_type));


 float_type Eo = real(Er[0]), En = 0; 
 for (int nt=1; nt < Nt; nt++)
 { 
  En = real(Er[nt]);
  for (int n=0; n<2*N; n++) { gamma[n*Nx] = 0.0; kappa[n*Nx] = 1.0; }
  for (int nx=1; nx<(Nx-1); nx++)
  {
   float_type x = (nx+1)*dx; 
   float_type A = dt/dx/dx - dt*(1.0/x - 1.0)/2.0/dx;
   float_type C = dt/dx/dx + dt*(1.0/x - 1.0)/2.0/dx; 
   for (int n=0; n<2*N; n++)
   {	
	float_type z = (n<N) ? Zeff[n] : (-Zeff[n-N]);
        float_type Enc = En/z/z/z;
	float_type Eoc = Eo/z/z/z; 
	f_complex B = f_complex(-2.0*dt/dx/dx + dt*x*Enc/2.0, 1.0); 
    	f_complex D = f_complex(-2.0*dt/dx/dx + dt*x*Eoc/2.0,-1.0); 
    
    	f_complex F = -A*f[nx-1 + n*Nx] - C*f[nx + 1 + n*Nx] - D*f[nx + n*Nx]; 
    
        f_complex G = (A*kappa[nx-1+n*Nx]+B);
        kappa[nx + n*Nx] = -C/G; 
        gamma[nx + n*Nx] = (F-A*gamma[nx -1 + n*Nx])/G;
   } 
  }
  
  for (int n=0; n<2*N; n++) f[Nx-1 + n*Nx] = gamma[Nx-2 + n*Nx]/((float_type)1.0-kappa[Nx-2 + n*Nx]);
  for (int nx=Nx-1; nx>0; nx--) for (int n=0; n<2*N; n++) f[nx-1 + n*Nx] = kappa[nx-1 + n*Nx]*f[nx + n*Nx] + gamma[nx-1 + n*Nx]; 
  Eo = En; 	

  if ((nt % tfac) == 0)
  {
   for (int i=0; i<3; i++) 
   {
 //   float_type x1 = dx, F1 = exp(-x1)*pow(x1,i);
    for (int n=0; n<2*N; n++) Ixn[i+3*n] = 0; 
    for (int nx=0; nx<Nx; nx++)
    {
     float_type x =   (nx+1)*dx, F1 = exp(-x)*pow(x, i);
     for (int n=0; n<2*N; n++) Ixn[i+3*n] += abs2(f[nx+n*Nx])*F1; 
     //F1 = F2; 
    }
    for (int n=0; n<2*N; n++) Ixn[i+3*n] *= dx; 
   }
   NL_SMALLBUFFER1[nt/tfac] = 0; 
   for (int n=0; n<N; n++) NL_SMALLBUFFER1[nt/tfac] += (Ixn[2+3*n]*Ixn[0+3*(n+N)] - Ixn[0+3*n]*Ixn[2+3*(n+N)])/(Ixn[1+3*n]*Ixn[0+3*(n+N)] + Ixn[0+3*n]*Ixn[1+3*(n+N)])*alpha[n]/Zeff[n]; 
   NL_SMALLBUFFER1[nt/tfac] *= 0.5*ELECTRON_CHARGE*BOHR_RADIUS*NEUTRAL_DENSITY;  
  }
 }
	
 free(f); free(gamma); free(kappa); free(Er); free(Ixn); 

#ifndef _UNWRAP_FREQUENCIES 
 NL_SMALLBUFFER1[0] = 0; for (int nt=1; nt<N_T; nt++) NL_SMALLBUFFER1[nt] *= exp(-j*OMEGA0*(TMIN+TSTEP*nt));
#endif 

 fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)output);
}


void calculate_plasmadensity_2float(f_complex* input, float_type* output, size_t N)
{
//#ifndef _ENABLE_CUDA
 for (size_t i=0; i<N; i++) calculate_plasmadensity_losses_small(input+N_T*i, output+N_T*i);
 return;
/*#else	

 f_complex* cuda_field_piece, *cuda_buf; float_type *cuda_plasma_piece;  

 size_t device_freemem=0, device_totalmem=0;
 int device_processN; MPI_Comm_size(DEVICE_COMM, &device_processN);

 (cudaMemGetInfo(&device_freemem, &device_totalmem));

 size_t cuda_pointN = exp2(floor(log2((double)device_freemem/N_T/sizeof(f_complex)/1.5/device_processN)));
 if (N < cuda_pointN) cuda_pointN = N;

 int    blocksize = 64;  if (blocksize > cuda_pointN) blocksize=cuda_pointN;
 size_t cuda_piece_size  = N_T*cuda_pointN*sizeof(f_complex);
 size_t cuda_piece_sizef = N_T*cuda_pointN*sizeof(float_type);

 (cudaMalloc((void**)&cuda_field_piece,       cuda_piece_size ));
 (cudaMalloc((void**)&cuda_buf,               cuda_piece_size ));
 (cudaMalloc((void**)&cuda_plasma_piece,      cuda_piece_sizef));

 for (long i=0; i<N; i+= cuda_pointN)
 {
  if ((N-i)<cuda_pointN) 
  {
	  cuda_pointN = N-i;
	  cuda_piece_size  = N_T*cuda_pointN*sizeof(f_complex);
	  cuda_piece_sizef = N_T*cuda_pointN*sizeof(float_type);
  }

  (cudaMemcpy(cuda_field_piece,    input+i*N_T,   cuda_piece_size, cudaMemcpyHostToDevice)) ;
  int blockN = cuda_pointN/blocksize, r = cuda_pointN % blocksize;
  calculate_plasma_2float_kernel<<<blockN,blocksize>>>(cuda_field_piece, cuda_plasma_piece, cuda_buf);
  if (r != 0) calculate_plasma_2float_kernel<<<1,r>>>(cuda_field_piece+N_T*blockN*blocksize, cuda_plasma_piece+N_T*blockN*blocksize, cuda_buf+N_T*blockN*blocksize);
  (cudaMemcpy(output+i*N_T, cuda_plasma_piece,    cuda_piece_sizef, cudaMemcpyDeviceToHost)) ;
 }
  
 (cudaFree(cuda_field_piece)); 
 (cudaFree(cuda_buf));
 (cudaFree(cuda_plasma_piece)); 
#endif*/
}

void calculate_maxplasmadensity_2float(f_complex* input, float_type* output, size_t N)
{
//#ifndef _ENABLE_CUDA

  float_type* b = (float_type*)NL_SMALLBUFFER1;
  for (size_t i=0; i<N; i++) 
  {
	  float_type ro=0;
	  calculate_plasmadensity_losses_small(input+N_T*i, b); 
	  for (size_t nt=0; nt<N_T; nt++) ro=max(ro, b[nt]);
	  output[i]=ro;
  }

/*#else
 f_complex* cuda_field_piece, *cuda_buf; float_type *cuda_maxplasma_piece;  

 size_t device_freemem=0, device_totalmem=0;
 int device_processN; MPI_Comm_size(DEVICE_COMM, &device_processN);

 (cudaMemGetInfo(&device_freemem, &device_totalmem));

 size_t cuda_pointN = exp2(floor(log2((double)device_freemem/(N_T+0.5)/sizeof(f_complex)/device_processN))); 
 if (cuda_pointN > N) cuda_pointN=N;

 int    blocksize = 64;  if (blocksize > cuda_pointN) blocksize=cuda_pointN;
 size_t cuda_piece_size  = N_T*cuda_pointN*sizeof(f_complex);
 size_t cuda_piece_sizef =     cuda_pointN*sizeof(float_type);

 (cudaMalloc((void**)&cuda_field_piece,       cuda_piece_size ));
 (cudaMalloc((void**)&cuda_buf,               cuda_piece_size ));
 (cudaMalloc((void**)&cuda_maxplasma_piece,   cuda_piece_sizef));

 for (long i=0; i<N; i+= cuda_pointN)
 {
  if ((N-i)<cuda_pointN) 
  {
	  cuda_pointN = N-i;
	  cuda_piece_size  = N_T*cuda_pointN*sizeof(f_complex);
	  cuda_piece_sizef =     cuda_pointN*sizeof(float_type);
  }

  (cudaMemcpy(cuda_field_piece,    input+i*N_T,   cuda_piece_size, cudaMemcpyHostToDevice)) ;
  int blockN = cuda_pointN/blocksize, r = cuda_pointN % blocksize;
  calculate_maxplasma_kernel<<<blockN,blocksize>>>(cuda_field_piece, cuda_maxplasma_piece, cuda_buf);
  if (r != 0) calculate_maxplasma_kernel<<<1,r>>>(cuda_field_piece+N_T*blockN*blocksize, cuda_maxplasma_piece+blockN*blocksize, cuda_buf+N_T*blockN*blocksize);
  (cudaMemcpy(output+i*N_T, cuda_maxplasma_piece,    cuda_piece_sizef, cudaMemcpyDeviceToHost)) ;
 }
  
 (cudaFree(cuda_field_piece)); 
 (cudaFree(cuda_maxplasma_piece)); 
 (cudaFree(cuda_buf));	
#endif*/
}





//----------------------------------------------
//------------------------------------------------


#ifdef _ENABLE_CUDA

void cuda_load_const()
{
 int deviceN = 0; 
 (cudaGetDeviceCount(&deviceN)); 
 if (deviceN == 0) throw "No CUDA-enabled devices found.";
 int devID_start = load_namedint(SC_FID, "DEVICE_ID_START", true, 0);
 int mydevice = 0;

#ifndef _GPU_CLUSTER_MAPPING
 if (ISMASTER) printf("\n %d CUDA devices found; process number = %d", deviceN, PROCESS_N); 
 
 if (PROCESS_N <= deviceN)  mydevice = (PROCESS_RANK+devID_start) % deviceN;
 else mydevice = (PROCESS_RANK*deviceN/PROCESS_N + devID_start)%deviceN;
#else
 mydevice = PROCESS_RANK%deviceN;
 if (ISMASTER) printf("\n %d CUDA devices per node", deviceN);
#endif

 (cudaSetDevice(mydevice)); 
 MPI_Comm_split(MPI_COMM_WORLD, mydevice, PROCESS_RANK, &DEVICE_COMM);
 int localrank = 0; MPI_Comm_rank(DEVICE_COMM, &localrank);

 printf("\n[%d]:Assigning to device %d, local rank %d.", PROCESS_RANK, mydevice, localrank);


 float_type* plasma_func_device, *raman_func_device;
 float_type* omega_device, *wavenum_device, *hodisp_device;

 float_type pf_re = real(PLASMA_FACTOR); 
 float_type pf_im = imag(PLASMA_FACTOR);

 (cudaMalloc(&plasma_func_device, sizeof(f_complex)*N_T)); 
 (cudaMalloc(&raman_func_device,  sizeof(f_complex)*N_T));

 (cudaMalloc(&wavenum_device,    sizeof(f_complex)*N_T));
 (cudaMalloc(&hodisp_device,     sizeof(float_type)*N_T));
 (cudaMalloc(&omega_device,      sizeof(float_type)*N_T));

 (cudaMemcpy(plasma_func_device, PLASMA_FUNC,        sizeof(f_complex)*N_T, cudaMemcpyHostToDevice));
 (cudaMemcpy(raman_func_device,  RAMAN_FUNCTION,     sizeof(f_complex)*N_T, cudaMemcpyHostToDevice));

 (cudaMemcpy(wavenum_device,  WAVENUMBER,     sizeof(f_complex) *N_T, cudaMemcpyHostToDevice));
 (cudaMemcpy(hodisp_device,   HO_DISPERSION,  sizeof(float_type)*N_T, cudaMemcpyHostToDevice));
 (cudaMemcpy(omega_device,    OMEGA,          sizeof(float_type)*N_T, cudaMemcpyHostToDevice));

 (cudaMemcpyToSymbol("cuda_plasma_factor_re",(void*)&pf_re,      sizeof(float_type),0,cudaMemcpyHostToDevice));
 (cudaMemcpyToSymbol("cuda_plasma_factor_im",(void*)&pf_im,      sizeof(float_type),0,cudaMemcpyHostToDevice));

 (cudaMemcpyToSymbol("cuda_plasma_func",  (void*)&plasma_func_device, sizeof(float_type*),0,cudaMemcpyHostToDevice));
 (cudaMemcpyToSymbol("cuda_raman_func",   (void*)&raman_func_device,  sizeof(float_type*),0,cudaMemcpyHostToDevice));
 (cudaMemcpyToSymbol("cuda_wavenum",      (void*)&wavenum_device,    sizeof(float_type*),0,cudaMemcpyHostToDevice));
 (cudaMemcpyToSymbol("cuda_omega",        (void*)&omega_device,      sizeof(float_type*),0,cudaMemcpyHostToDevice));
 (cudaMemcpyToSymbol("cuda_hodisp",        (void*)&hodisp_device,    sizeof(float_type*),0,cudaMemcpyHostToDevice));

  (cudaMemcpyToSymbol("cuda_N_T",  (void*)&N_T,  sizeof(int),       0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_TMIN", (void*)&TMIN,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_TMAX", (void*)&TMAX,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  float_type tstep = (TMAX-TMIN)/N_T; 
  (cudaMemcpyToSymbol("cuda_TSTEP",(void*)&tstep, sizeof(float_type),0,cudaMemcpyHostToDevice)); 

  (cudaMemcpyToSymbol("cuda_IONIZATION_POTENTIAL",  (void*)&IONIZATION_POTENTIAL,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_AVALANCHE_CROSSSECTION", (void*)&AVALANCHE_CROSSSECTION,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_RECOMBINATION_TAU", (void*)&RECOMBINATION_TAU,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_GROUP_VELOCITY", (void*)&GROUP_VELOCITY,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_RECOMBINATION_TAU", (void*)&RECOMBINATION_TAU,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_AMBIENT_CARRIER_DENSITY", (void*)&AMBIENT_CARRIER_DENSITY,  sizeof(float_type),0,cudaMemcpyHostToDevice));

  (cudaMemcpyToSymbol("cuda_N2", (void*)&NONLIN_REFRINDEX ,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_N4", (void*)&NONLIN_REFRINDEX4,  sizeof(float_type),0,cudaMemcpyHostToDevice));

  (cudaMemcpyToSymbol("cuda_PONDEROMOTIVE_COEFFICIENT", (void*)&PONDEROMOTIVE_COEFFICIENT,  sizeof(float_type),0,cudaMemcpyHostToDevice));

#ifdef MULTI_LEVEL_IONIZATION 
  (cudaMemcpyToSymbol("cuda_IONIZATION_LEVEL_N",     (void*)&IONIZATION_LEVEL_N, sizeof(int), 0, cudaMemcpyHostToDevice)); 
  
  float_type* Us_device;  (cudaMalloc(&Us_device,    IONIZATION_LEVEL_N*sizeof(float_type))); 
  (cudaMemcpyToSymbol("cuda_IONIZATION_POTENTIALS",   (void*)&Us_device, sizeof(float_type*), 0, cudaMemcpyHostToDevice)); 
  (cudaMemcpy(Us_device, IONIZATION_POTENTIALS,      IONIZATION_LEVEL_N*sizeof(float_type), cudaMemcpyHostToDevice)); 

  if (N_T < IONIZATION_LEVEL_N*2) throw "cuda_load_const(): N_T should be greater than 2 IONIZATION_LEVEL_N, because one of nonlinear-response buffers is used to store ion densities and ionization rates during plasma density calculation.";
  int ionization_levels_n = IONIZATION_LEVEL_N;
#else 
  int ionization_levels_n  = 1;  
#endif
  float_type* Wln_device; (cudaMalloc(&Wln_device, IONIZATION_N*ionization_levels_n*sizeof(float_type))); 
  (cudaMemcpyToSymbol("cuda_IONIZATION_RATE_LN", (void*)&Wln_device, sizeof(float_type*),0, cudaMemcpyHostToDevice));
  (cudaMemcpy(Wln_device, IONIZATION_RATE_LN, IONIZATION_N*ionization_levels_n*sizeof(float_type), cudaMemcpyHostToDevice));

  float_type kxstep = 2*M_PI/(XMAX-XMIN), kystep = 2*M_PI/(YMAX-YMIN);
 
  (cudaMemcpyToSymbol("cuda_kxstep", (void*)&kxstep,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_kystep", (void*)&kystep,  sizeof(float_type),0,cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_OMEGA_MAX", (void*)&OMEGA_MAX, sizeof(float_type), 0, cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_OMEGA_MIN", (void*)&OMEGA_MIN, sizeof(float_type), 0, cudaMemcpyHostToDevice));
 
  (cudaMemcpyToSymbol("cuda_NEUTRAL_DENSITY", (void*)&NEUTRAL_DENSITY, sizeof(float_type), 0, cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_RAMAN_FRACTION", (void*)&RAMAN_FRACTION, sizeof(float_type), 0, cudaMemcpyHostToDevice));

  (cudaMemcpyToSymbol("cuda_wavenum0", (void*)&WAVENUMBER0, sizeof(float_type), 0, cudaMemcpyHostToDevice));
  (cudaMemcpyToSymbol("cuda_omega0",   (void*)&OMEGA0,      sizeof(float_type), 0, cudaMemcpyHostToDevice));

#ifdef MULTIPHOTON_IONIZATION
 (cudaMemcpyToSymbol("cuda_K_MPI", (void*)&K_MPI, sizeof(int), 0, cudaMemcpyHostToDevice));
 (cudaMemcpyToSymbol("cuda_BETA_MPI_LN", (void*)&BETA_MPI_LN, sizeof(float_type), 0, cudaMemcpyHostToDevice));
#endif

#ifdef THIRD_HARMONICS
 (cudaMemcpyToSymbol("cuda_th_factor", (void*)&TH_FACTOR, sizeof(float_type), 0, cudaMemcpyHostToDevice)); 
#endif
}


void cuda_free_const()
{
 float_type* pplasma_func, **praman_func, *pwavenum, *pomega, *phodisp;

 (cudaMemcpyFromSymbol(&pplasma_func, "cuda_plasma_func",  sizeof(float_type*), 0, cudaMemcpyDeviceToHost)); cudaFree(pplasma_func);
 (cudaMemcpyFromSymbol(&praman_func,  "cuda_raman_func",   sizeof(float_type*), 0, cudaMemcpyDeviceToHost)); cudaFree(praman_func);
 (cudaMemcpyFromSymbol(&pwavenum,     "cuda_wavenum",      sizeof(float_type*), 0, cudaMemcpyDeviceToHost)); cudaFree(pwavenum);
 (cudaMemcpyFromSymbol(&pomega,       "cuda_omega",        sizeof(float_type*), 0, cudaMemcpyDeviceToHost)); cudaFree(pomega);
 (cudaMemcpyFromSymbol(&phodisp,      "cuda_hodisp",       sizeof(float_type*), 0, cudaMemcpyDeviceToHost)); cudaFree(phodisp);

 if (cuda_use_pinned_memory)
 {
  cudaHostUnregister(FIELD); cudaHostUnregister(BIGBUFFER1); cudaHostUnregister(BIGBUFFER2);
 }

 if (!cuda_do_extra_memtransfer)
 {
  cudaFree(cuda_field); cudaFree(cuda_buf1); cudaFree(cuda_buf2);
  cudaFree(cuda_bufs);
 }
 
 
}


__global__ void calculate_Lresponse_kernel(float_type* field, float_type* buf1, float_type* kt2_, float_type ZSTEP, size_t N_T, size_t n)
{
  size_t  nlocal = blockDim.x*blockIdx.x+threadIdx.x;
 

  for (size_t nw=0; nw<N_T; nw++)
  {
    size_t ofs = 2*(nw + N_T*nlocal); 
#ifndef NO_SPACE_TIME_FOCUSING
     float_type k1  = cuda_wavenum[2*nw];
#else
     float_type k1 = cuda_wavenum0;
#endif	
     float_type k2 = cuda_wavenum[2*nw+1];
     float_type w = cuda_omega[nw]; 
	 float_type kt2n = kt2_[nlocal+n]/k1/k1;
#ifndef NO_DIFFRACTION
    if (w < 0 || kt2n > MAX_KT2) { for (size_t i=0; i<2; i++) field[ofs+i]=0.0; continue;}
#else
    if (w < 0) { for (size_t i=0; i<2; i++) field[ofs+i]=0.0; continue;}
#endif
     
     float_type reE = field[ofs], imE = field[ofs+1];
#ifdef PARABOLICAL_DIFRACTION

	
#ifndef NO_DIFFRACTION
     	float_type phi = ZSTEP*(k1*kt2n/(float_type)2.0 - cuda_hodisp[nw]);
#else
	 	float_type phi = ZSTEP*(-cuda_hodisp[nw]);
#endif
	 float_type M   = exp_f(ZSTEP*k2);
	 float_type c,s;
	 sincos_p(phi,&s,&c);

#ifdef NONLINEARITY_ON      
     reE += ZSTEP*buf1[ofs]; 
     imE += ZSTEP*buf1[ofs+1];
     field[ofs]   = (reE*c - imE*s)*M;// + ZSTEP*buf1[ofs];
     field[ofs+1] = (imE*c + reE*s)*M;// + ZSTEP*buf1[ofs+1];
#endif

#else
	 float_type sqHO = device_sqrtHO(-kt2n); 
	 float_type phi = ZSTEP*(-cuda_hodisp[nw] - k1*sqHO);
	 float_type M   = exp_f(ZSTEP*k2);  
	 float_type c,s;
	 sincos_p(phi,&s,&c);
#ifdef NONLINEARITY_ON      
     reE += ZSTEP*buf1[ofs]  /(1.0+sqHO);   
     imE += ZSTEP*buf1[ofs+1]/(1.0+sqHO); 
     field[ofs]   = (reE*c - imE*s)*M;// + ZSTEP*buf1[ofs]  /(1+sqHO);
     field[ofs+1] = (imE*c + reE*s)*M;// + ZSTEP*buf1[ofs+1]/(1+sqHO);
#endif

#endif


#ifndef NONLINEARITY_ON
     field[ofs]   = (reE*c - imE*s)*M;
     field[ofs+1] = (imE*c + reE*s)*M;
#endif
   }
}


__global__ void calculate_NLresponse_kernel_strided(f_complex* input, f_complex* output, float_type zstep, int N_T, float_type* bufs, float_type* maxA2buf, float_type* maxNL2buf)
{
  
 int ofs =          N_T*(blockDim.x*blockIdx.x   + threadIdx.x);
#ifdef _ODE_RK4
 int bufofs =      20*N_T*(blockDim.x*blockIdx.x)  + threadIdx.x; 
#endif
#ifdef _ODE_HEUN
 int bufofs =      16*N_T*(blockDim.x*blockIdx.x)  + threadIdx.x; 
#endif
 #ifdef _ODE_EULER
 int bufofs =      14*N_T*(blockDim.x*blockIdx.x)  + threadIdx.x; 
#endif


 float_type* inp  = (float_type*)(input +ofs);
 float_type* outp = (float_type*)(output+ofs);


 int stride = blockDim.x;
 int stridei = 31-__clz(stride);
 int ntstride = N_T * 2*stride ;

 float_type* b   = bufs+bufofs;
 float_type* input_buf   = b+2*N_T*stride;
 float_type* output_buf1 = b+4*N_T*stride;

#ifdef _ODE_RK4
 float_type* temp_bufs   = b+12*N_T*stride;
#endif
#ifdef _ODE_HEUN
 float_type* temp_bufs   = b+8*N_T*stride;
#endif
#ifdef _ODE_EULER
 float_type* temp_bufs   = b+6*N_T*stride;
#endif


 float_type maxI_ = 0, maxNL2_ = 0;
 //Runge-Kutta 4-step:
 for (unsigned int i=0; i<2*N_T; i++) 
 {
	 b[i << stridei] = inp[i];
 }

 for (int i=0; i<2*N_T; i++) input_buf[i*stride]=b[i*stride];

 for (int i=0; i<N_T; i++) maxI_ = max_p(maxI_, cuda_abs2(b[ 2*i   *stride], b[(2*i+1)*stride]));

#ifdef _ODE_RK4
 float_type runge_kutta_factor[3] = {0.5,0.5,1};
 for (int nr=0; nr<4; nr++)
 {
    if (nr > 0) for (int i=0; i<2*N_T; i++) input_buf[i*stride] = b[i*stride] + runge_kutta_factor[nr-1]*zstep*output_buf1[i*stride+(nr-1)*ntstride];
	calculate_single_NLresponse_kernel_strided(input_buf, output_buf1+nr*ntstride, N_T, temp_bufs);
 }
#endif 

#ifdef _ODE_HEUN
 for (int nr=0; nr<2; nr++)
 {
    if (nr > 0) for (int i=0; i<2*N_T; i++) input_buf[i*stride] = b[i*stride] + zstep*output_buf1[i*stride+(nr-1)*ntstride];
	calculate_single_NLresponse_kernel_strided(input_buf, output_buf1+nr*ntstride, N_T, temp_bufs);
 }
#endif 

#ifdef _ODE_EULER
 calculate_single_NLresponse_kernel_strided(input_buf, output_buf1, N_T, temp_bufs);
#endif 

 for (unsigned int nt=0; nt<N_T; nt++)
 {
#ifdef _ODE_EULER
  float_type reout = output_buf1[2*nt*stride]; float_type imout=output_buf1[(2*nt+1)*stride]; 
#endif
#ifdef _ODE_HEUN
  float_type reout = (output_buf1[ 2*nt   *stride] + output_buf1[ 2*nt   *stride + ntstride])/2.0; 
  float_type imout = (output_buf1[(2*nt+1)*stride] + output_buf1[(2*nt+1)*stride + ntstride])/2.0; 
#endif
#ifdef _ODE_RK4
  float_type reout = (output_buf1[ 2*nt   *stride] + 2.0*output_buf1[ 2*nt   *stride + ntstride] + 2.0*output_buf1[ 2*nt   *stride + 2*ntstride] + output_buf1[ 2*nt   *stride + 3*ntstride])/6.0; 
  float_type imout = (output_buf1[(2*nt+1)*stride] + 2.0*output_buf1[(2*nt+1)*stride + ntstride] + 2.0*output_buf1[(2*nt+1)*stride + 2*ntstride] + output_buf1[(2*nt+1)*stride + 3*ntstride])/6.0; 
#endif
  maxNL2_ = max_p(maxNL2_, cuda_abs2(reout, imout)); 
  outp[2*nt] = reout; outp[2*nt+1]=imout; 
 }

 maxA2buf[blockIdx.x*blockDim.x + threadIdx.x] = maxI_; 
 maxNL2buf[blockIdx.x*blockDim.x + threadIdx.x] = 2.31; //maxNL2_; 
  
}



__device__ void calculate_single_NLresponse_kernel_strided(float_type* input, float_type* output, int N_T, float_type* bufs)
{
 int stride = blockDim.x;
#ifdef PLASMA_FULL_DISPERSION
 for (int nw=0; nw<2*N_T; nw++) bufs[(2*N_T+nw)*stride] = input[nw*stride];
#endif
 fft_device_strided(input, N_T, 1, stride);
 calculate_plasmadensity_losses_small_device(input, bufs, stride, 2, bufs+4*N_T*stride, bufs+6*N_T*stride); 

 
#ifndef NO_PLASMARESPONSE
#ifndef PLASMA_DISPERSION
#ifndef PLASMA_FULL_DISPERSION 
for (int nt=0; nt<N_T; nt++)
 {
   float_type reE = input[2*nt*stride], imE = input[(2*nt+1)*stride];
   float_type ro  = bufs[(2*nt)*stride];

   float_type re_Fro  = -cuda_plasma_factor_re*ro;
   float_type im_Fro  = -cuda_plasma_factor_im*ro;
   float_type resFro, imsFro;
   if (fabs(re_Fro) < 0.1) device_sqrtHO(re_Fro,im_Fro, &resFro, &imsFro);
   else {sqrtc(1+re_Fro, im_Fro, &resFro, &imsFro);  resFro -= 1;}

   float_type reR  =  cuda_wavenum0* imsFro;
   float_type imR  = -cuda_wavenum0* resFro; 
 
   bufs[(4*N_T+2*nt)  *stride]  += reR*reE - imR*imE;
   bufs[(4*N_T+2*nt+1)*stride]  += reR*imE + imR*reE;
 }
#else
 for (int nt=0; nt<N_T; nt++)
 {
   float_type ro  = bufs[(2*nt)*stride]; 
   float_type reM = (float_type)1.0/(float_type)N_T;
   float_type imM = 0; 
   float_type rew, imw;  sincos_f((float_type)2.0*M_PI*(float_type)nt/(float_type)N_T, &imw, &rew); 
   bufs[(4*N_T+2*nt)*stride] = 0; bufs[(4*N_T+2*nt+1)*stride]=0;

   for (int nw=0; nw<N_T; nw++)
   {
   	float_type reT1, imT1, reT2, imT2;
	reT1  = cuda_plasma_func[2*nw]*ro; imT1  = cuda_plasma_func[2*nw+1]*ro;
 	if (fabs(reT1) < 0.1) device_sqrtHO(reT1,imT1, &reT2, &imT2);
   	else {sqrtc(1+reT1, imT1, &reT2, &imT2);  reT2 -= 1.0;} 		
	reT1 = cuda_wavenum[2*nw]; imT1 = reT1*imT2; reT1 = reT1*reT2;    //T1 = k[w];  T1 = T1*T2 => T1 = T2*k(w

	float_type reAw = bufs[(2*N_T+2*nw)*stride], imAw = bufs[(2*N_T+2*nw+1)*stride];  

	reT2 = reT1*reAw - imT1*imAw; imT2 = reT1*imAw + imT1*reAw; //T2 = T1*A(w)  => T2 = T2*k(w)*A(w)

  	
	bufs[(4*N_T+2*nt)*stride]   +=  (reT2*imM + imT2*reM);      //bufs += -i*k(w)*T2*A(w)*M;  where M = exp(i*2*pi*nt*nw/N_T); 
	bufs[(4*N_T+2*nt+1)*stride] += -(reT2*reM - imT2*imM);

	reT1 = reM; imT1 = imM;		  // M=M*exp(i*2*pi*nt/N_T); 
	reM = reT1*rew - imT1*imw; 
	imM = reT1*imw + imT1*rew; 
   }	
 }
#endif

 fft_device_strided(bufs+4*N_T*stride, N_T, -1, stride);


#else
  for (int nt=0; nt<N_T; nt++)
  {
   float_type reE = input[2*nt*stride], imE = input[(2*nt+1)*stride];
   float_type ro  = bufs[(2*nt)*stride];
   
   bufs[(2*N_T+2*nt)  *stride] = ro*reE;
   bufs[(2*N_T+2*nt+1)*stride] = ro*imE;
  }

 

  fft_device_strided(bufs+4*N_T*stride, N_T, -1,stride); fft_device_strided(bufs+2*N_T*stride, N_T, -1, stride);

  for (int nw=0; nw<N_T; nw++) 
  {
   bufs[(4*N_T+2*nw)  *stride] += cuda_plasma_func[2*nw]*bufs[(2*N_T+2*nw)  *stride] - cuda_plasma_func[2*nw+1]*bufs[(2*N_T+2*nw+1)*stride];
   bufs[(4*N_T+2*nw+1)*stride] += cuda_plasma_func[2*nw]*bufs[(2*N_T+2*nw+1)*stride] + cuda_plasma_func[2*nw+1]*bufs[(2*N_T+2*nw)  *stride];
  }
#endif

#else
 for (int nw=0; nw<2*N_T; nw++) bufs[(4*N_T+nw)  *stride] = 0; 
#endif

 for (int nt=0; nt<N_T; nt++) {bufs[ 2*nt   *stride] = cuda_abs2(input[2*nt*stride], input[(2*nt+1)*stride]); bufs[(2*nt+1)*stride] = 0;}
 fft_device_strided(bufs, N_T, -1, stride);

 for (int nw=0; nw<N_T; nw++)
 {
  float_type reIw = bufs[2*nw*stride], imIw = bufs[(2*nw+1)*stride];
  bufs[(2*N_T+2*nw)  *stride] = cuda_raman_func[2*nw]*reIw - cuda_raman_func[2*nw+1]*imIw; 
  bufs[(2*N_T+2*nw+1)*stride] = cuda_raman_func[2*nw]*imIw + cuda_raman_func[2*nw+1]*reIw;
 }
 fft_device_strided(bufs+2*N_T*stride, N_T,  1, stride); 

 for (int nt=0; nt<N_T; nt++)
 {
  float_type reE  = input[2*nt*stride], imE = input[(2*nt+1)*stride]; 
  float_type I    = cuda_abs2(reE, imE); 
#ifdef THIRD_HARMONICS 
#ifndef _UNWRAP_FREQUENCIES
  float_type carrier_phase = 2.0*cuda_omega0*(cuda_TMIN + nt*cuda_TSTEP);
#else  
  float_type carrier_phase = 0; 
#endif 

  float_type s,c;     sincos_p(carrier_phase,&s,&c);
  float_type s2, c2;  sincos_p(2.0*carrier_phase, &s2, &c2);
  float_type reE2 = reE*reE - imE*imE; 
  float_type imE2 = 2.0*reE*imE;
  float_type reE4 = reE2*reE2 - imE2*imE2;
  float_type imE4 = 2.0*reE2*imE2; 

  float_type reRf = cuda_N2*(bufs[(2*N_T+2*nt)*stride]   + (1-cuda_RAMAN_FRACTION)*cuda_th_factor/3.0*(reE2*c - imE2*s)) + cuda_N4*(I*I + I*(reE2*c - imE2*s)/2.0 + (reE4*c2 - imE4*s2)/10.0);
  float_type imRf = cuda_N2*(bufs[(2*N_T+2*nt+1)*stride] + (1-cuda_RAMAN_FRACTION)*cuda_th_factor/3.0*(reE2*s + imE2*c)) + cuda_N4*(      I*(reE2*s + imE2*c)/2.0 + (reE4*s2 - imE4*c2)/10.0);
#else 
  float_type reRf = cuda_N2*bufs[(2*N_T+2*nt)*stride] + cuda_N4*I*I; 
  float_type imRf = cuda_N2*bufs[(2*N_T+2*nt+1)*stride] + cuda_N4*I*I;
#endif 

  bufs[ 2*nt   *stride] = reRf*reE - imRf*imE; 
  bufs[(2*nt+1)*stride] = imRf*reE + reRf*imE;
 }
 fft_device_strided(bufs, N_T, -1, stride);

 for (int nw=0; nw<N_T; nw++)
 { 
#ifndef NO_SHOCK
  float_type w = cuda_omega[nw];
#else
  float_type w = cuda_omega0;
#endif
  float_type reRf = bufs[2*nw*stride], imRf = bufs[(2*nw+1)*stride]; 
  bufs[(2*N_T+2*nw)  *stride]  =   imRf*w/LIGHT_VELOCITY;
  bufs[(2*N_T+2*nw+1)*stride]  =  -reRf*w/LIGHT_VELOCITY;
 }


#ifdef NONLINEARITY_ON
 for (int i=0; i<2*N_T; i++)
 {
   float_type w = cuda_omega[i>>1]; if (w<0) {output[i*stride]=0.0; continue;} 
   output[i*stride]=(bufs[(4*N_T+i)*stride]+bufs[(2*N_T+i)*stride])*device_hfgaussfilter(w)*device_lfgaussfilter(w);
 }
#else
 for (int i=0; i<2*N_T; i++) output[i*stride]=0.0;
#endif
//for (int i=0; i<2*N_T; i++) output[i*stride]=bufs[(2*N_T+i)*stride];
}



__device__ inline float_type plasma_source_function_device(float_type reA, float_type imA, float_type ro)
{
	return photoionization_function_device(reA, imA, ro)  + avalanche_ionization_function_device(reA, imA, ro) - recombination_function_device(ro);
}


__device__ inline float_type avalanche_ionization_function_device(float_type reA, float_type imA, float_type ro)
{
	if (ro < 0 || ro > cuda_NEUTRAL_DENSITY) return 0;
	float_type I = cuda_abs2(reA, imA); 
    return cuda_AVALANCHE_CROSSSECTION/(cuda_IONIZATION_POTENTIAL+cuda_PONDEROMOTIVE_COEFFICIENT*I)*I*ro*(1-ro/cuda_NEUTRAL_DENSITY);
}

__device__ inline float_type recombination_function_device(float_type ro)
{
	return (ro>0)?(ro/cuda_RECOMBINATION_TAU):0;
}

__device__ inline float_type photoabsorbtion_function_device(float_type reA, float_type imA, float_type ro)
{
	float_type I = cuda_abs2(reA,imA);
	return ((I > IONIZATION_MIN_I/INTENSITY_DENOM)?((float_type)0.5*photoionization_function_device(reA, imA,ro)*((cuda_IONIZATION_POTENTIAL+cuda_PONDEROMOTIVE_COEFFICIENT*I)/I)):(0));
}

__device__ inline void calculate_plasmadensity_small_device(float_type* field, float_type* pro, float_type* buf)
{
        calculate_plasmadensity_small_device_strided(field, pro, 1, buf);
}

__device__ inline void calculate_plasmadensity_small_device_strided(float_type* field, float_type* pro, int stride, float_type* buf)
{
   calculate_plasmadensity_losses_small_device(field, pro, stride, 2, NULL, buf); 
}

__device__  void calculate_plasmadensity_losses_small_device(float_type* field, float_type* pro, int stride, int rostride, float_type* loss, float_type* buf)
{
    float_type tstep = cuda_TSTEP*IONRATE_DENOM;
#ifndef MULTI_LEVEL_IONIZATION
    if (pro) pro[0]=0; 
	if (loss) {loss[0]=0; loss[stride]=0;}
    float_type ro1 = 0;
	float_type ro2 = 0;

	float_type W10, W11, W20, W21, k1, k2;
#else
    float_type* W = buf + cuda_IONIZATION_LEVEL_N*stride; 
    float_type* ro = buf;  for (int i=0; i<cuda_IONIZATION_LEVEL_N; i++) ro[i*stride]=0;
    float_type  tro = 0, tl = 0;
    if (pro) pro[0] = 0; 
    if (loss) {loss[0]=0; loss[stride]=0;}
#endif 

    float_type reE0 = field[0],        imE0 = field[stride];
    float_type reE1=0, imE1=0, I0=cuda_abs2(reE0, imE0), I1=0;
      
   for (int nt=1; nt<cuda_N_T;   nt++)
   {
	  // Plasma density calculation in performed via Runge-Khutta method of second order (Heun method)
	  reE1 = field[(2*nt  )*stride];
	  imE1 = field[(2*nt+1)*stride];
	  I1 = cuda_abs2(reE1, imE1);
//	  printf("\nE1=%g+i%g, nt=%d, stride=%d, I1 = %g",reE1, imE1, nt, stride, I1);
#ifndef MULTI_LEVEL_IONIZATION
	  W10 = photoionization_function_device(reE0, imE0, ro1);    k1 = (tstep)*(W10 + avalanche_ionization_function_device(reE0, imE0, ro1) - recombination_function_device(ro1));	    
	  W11 = photoionization_function_device(reE1, imE1, ro1+k1); k2 = (tstep)*(W11 + avalanche_ionization_function_device(reE1, imE1, ro1+k1) - recombination_function_device(ro1+k1));	    
	  ro1 += 0.5*(k1+k2); if (ro1 > cuda_NEUTRAL_DENSITY) ro1 = cuda_NEUTRAL_DENSITY;

	  
	  if (pro)  pro[rostride*nt*stride] = ro1;
	  k1 = (I1>IONIZATION_MIN_I)?(-0.5*(W10*(cuda_IONIZATION_POTENTIAL+cuda_PONDEROMOTIVE_COEFFICIENT*I0)/I0)):0;
 
	  if (loss) {loss[(2*nt-2)*stride] = reE1*k1; loss[(2*nt-1)*stride]=imE1*k1;}
#else
	 if (loss) tl=0; 
	 float_type k1=0, k2=0;
	 if (I0 > IONIZATION_MIN_I/INTENSITY_DENOM)
	 {

	   photoionization_functionsN_device(reE0, imE0,     W,                                stride);
	   photoionization_functionsN_device(reE1, imE1,     W+cuda_IONIZATION_LEVEL_N*stride, stride);
	   k1 = tstep*W[0]*(cuda_NEUTRAL_DENSITY-ro[0]); k2 = tstep*W[cuda_IONIZATION_LEVEL_N*stride]*(cuda_NEUTRAL_DENSITY-ro[0]-k1); 
	   if (k1 < 0) k1=0; if (k2 < 0) k2 = 0;
	   ro[0] += 0.5*(k1+k2); if (ro[0] > cuda_NEUTRAL_DENSITY) ro[0]=cuda_NEUTRAL_DENSITY; 
	   tro = ro[0]; 
	
           tl = -(float_type)0.5*k1/tstep*(cuda_IONIZATION_POTENTIALS[0]+cuda_PONDEROMOTIVE_COEFFICIENT*I0)/I0;
	 
	   for (int n=1; n<cuda_IONIZATION_LEVEL_N; n++) 
	   {
 	      k1 = tstep*W[n*stride]*(ro[(n-1)*stride]-ro[n*stride]); k2 = tstep*W[(n+cuda_IONIZATION_LEVEL_N)*stride]*(ro[(n-1)*stride]-ro[n*stride]-k1); 
	      if (k1 < 0) k1=0; if (k2 < 0) k2 = 0;
	      ro[n*stride] = ro[n*stride] + 0.5*(k1+k2); 
	      if (ro[n*stride]>ro[(n-1)*stride]) ro[n*stride]=ro[(n-1)*stride];
	      tro  += ro[n*stride];
              tl -= (float_type)0.5*k1/tstep*(cuda_IONIZATION_POTENTIALS[n]+cuda_PONDEROMOTIVE_COEFFICIENT*I0)/I0;
	   }
	  }
  	  if (pro) {pro[rostride*stride*nt]=tro;}
          if (loss) {loss[(2*nt-2)*stride] = tl*reE0; loss[(2*nt-1)*stride]=tl*imE0;}
#endif
	  reE0=reE1; imE0 = imE1;
	  I0 = I1;
   } 
   if (loss) {loss[(2*cuda_N_T-2)*stride]=0; loss[(2*cuda_N_T-1)*stride]=0;}
}

__device__ inline void calculate_plasmadensity_small_device_strided_2float(float_type* field, float_type* pro, int stride, float_type* buf)
{
   calculate_plasmadensity_losses_small_device(field, pro, stride, 1, NULL, buf);
}

__device__ inline float_type calculate_maxplasmadensity_small_device(float_type* field, float_type* buf)
{
 //calculate_plasmadensity_losses_small_device(field, NULL, blockDim.x*gridDim.x, 1, NULL, buf);
 return 0;//buf[cuda_N_T-1];
}


__device__ float_type photoionization_function_device(float_type reA, float_type imA, float_type ro)
{
        float_type I = cuda_abs2(reA, imA); 
		if (ro > cuda_NEUTRAL_DENSITY) return 0;
        if (I > (IONIZATION_MIN_I/INTENSITY_DENOM))
	{
	     float_type lnI = log_p(I)+INTENSITY_DENOM_LN;
#ifdef MULTIPHOTON_IONIZATION
         return 0;
	//return ((ro<cuda_NEUTRAL_DENSITY)?(exp_p(cuda_K_MPI*lnI+cuda_BETA_MPI_LN-IONRATE_DENOM_LN)*(cuda_NEUTRAL_DENSITY-ro)):0);
#else

         int n = (float_type)(floor((lnI - IONIZATION_MIN_I_LN)/IONIZATION_I_LN_TOLERANCE));
         if (n < 0 || n > (IONIZATION_N-2)) return 0; 
         float_type lnIn  = IONIZATION_MIN_I_LN+n*IONIZATION_I_LN_TOLERANCE;

         float_type lnWn  = cuda_IONIZATION_RATE_LN[n];
         float_type lnWnp = cuda_IONIZATION_RATE_LN[n+1];
#ifdef IONIZATION_LINEAR_INTERP
#ifdef IONIZATION_GAS
		 float_type Wn = exp_p(lnWn+log_p(cuda_NEUTRAL_DENSITY-ro)), Wnp = exp(lnWnp+log_p(cuda_NEUTRAL_DENSITY-ro));
#else
		 float_type Wn = exp_p(lnWn), Wnp = exp_p(lnWnp);
#endif
		 float_type In = exp_p(lnIn), Inp = exp_p(lnIn+IONIZATION_I_LN_TOLERANCE);
		 return Wn + (I-In)*(Wnp-Wn)/(Inp-In);
#else
         float_type lnW = lnWn + (lnI-lnIn)*(lnWnp-lnWn)/IONIZATION_I_LN_TOLERANCE;

#ifdef IONIZATION_GAS
		 return exp_p(lnW + log_p(cuda_NEUTRAL_DENSITY-ro));
#else
         return exp_p(lnW); 
#endif
#endif
	}
	return 0;
#endif
}


#ifdef MULTI_LEVEL_IONIZATION

__device__ void photoionization_functionsN_device(float_type reA, float_type imA, float_type* W, int stride)
{

        float_type I = cuda_abs2(reA, imA);
        if (I < (IONIZATION_MIN_I/INTENSITY_DENOM)) {for (int i=0; i<cuda_IONIZATION_LEVEL_N; i++) W[i*stride]=0; return;}

	float_type lnI = log(I)+INTENSITY_DENOM_LN;
         int n = (int)floor((lnI - IONIZATION_MIN_I_LN)/IONIZATION_I_LN_TOLERANCE);
        if (n > IONIZATION_N-2 || n < 0) {for (int i=0; i<cuda_IONIZATION_LEVEL_N; i++) W[i*stride]=0; return;}  
	
#ifdef YUDIN_IVANOV_CORRECTION
         float_type theta = atan(imA/reA);  
         float_type wl = cuda_omega0*PLANCK_CONSTANT_REDUCED/HARTREE_ENERGY; 
	 float_type phase_v = cuda_omega0/cuda_wavenum0;
         float_type E = exp(lnI/2)*sqrt(2*VACUUM_PERMEABILITY*phase_v*phase_v/cuda_GROUP_VELOCITY);
         float_type T  = (E*E/ATOMIC_FIELD/ATOMIC_FIELD)/wl/wl/wl;
#endif
       
	for (int i=0; i<cuda_IONIZATION_LEVEL_N; i++)
	{
         float_type lnIn  = IONIZATION_MIN_I_LN+n*IONIZATION_I_LN_TOLERANCE;
 	 float_type lnWn  = cuda_IONIZATION_RATE_LN[IONIZATION_N*(i)+n];
         float_type lnWnp = cuda_IONIZATION_RATE_LN[IONIZATION_N*(i)+n+1];
	 float_type lnW = lnWn + (lnI-lnIn)*(lnWnp-lnWn)/IONIZATION_I_LN_TOLERANCE;

	 float_type F = 1.0; 
  #ifdef PPT_IONIZATION
  #ifdef YUDIN_IVANOV_CORRECTION
	float_type g = sqrt(2.0*cuda_IONIZATION_POTENTIALS[i]*INTENSITY_DENOM/IONRATE_DENOM/ELECTRON_CHARGE/ELECTRON_CHARGE*ELECTRON_MASS)*cuda_omega0/E;

	if (g < MAX_YI_GAMMA) F = exp(-T*(YI_Phi_device(theta, g)-YI_Phi_device(0,g)));         
  #endif
  #endif

	 W[i*stride] = exp(lnW)*F; 
	}
}
#endif


__device__  float_type photoionization_function_device2(float_type reA, float_type imA, float_type ro1, float_type ro2)
{
        float_type I = cuda_abs2(reA, imA); 
		if (ro2 > ro1) return 0;
		if (ro1 < 1e-5*cuda_NEUTRAL_DENSITY) return 0;
        if (I > (IONIZATION_MIN_I/INTENSITY_DENOM))
	{
	     float_type lnI = log_p(I)+INTENSITY_DENOM_LN;

         int n = (float_type)(floor((lnI - IONIZATION_MIN_I_LN)/IONIZATION_I_LN_TOLERANCE));
         if (n < 0 || n > (IONIZATION_N-2)) return 0; 
         float_type lnIn  = IONIZATION_MIN_I_LN+n*IONIZATION_I_LN_TOLERANCE;

         float_type lnWn  = cuda_IONIZATION_RATE_LN[IONIZATION_N+n];
         float_type lnWnp = cuda_IONIZATION_RATE_LN[IONIZATION_N+n+1];

		  float_type lnW = lnWn + (lnI-lnIn)*(lnWnp-lnWn)/IONIZATION_I_LN_TOLERANCE;

#ifdef IONIZATION_GAS
		 return exp_p(lnW + log_p(ro1-ro2));
#else
         return exp_p(lnW); 
#endif
	}
	return 0;
}




__global__ void calculate_plasma_2float_kernel(f_complex* input, float_type* output, f_complex* buf)
{
	int ofs = cuda_N_T*(blockDim.x*blockIdx.x   + threadIdx.x);
	calculate_plasmadensity_small_device_strided_2float((float_type*)(input+ofs), output+ofs, 1, (float_type*)(buf+ofs));
}


__global__ void calculate_maxplasma_kernel   (f_complex* input, float_type* output, f_complex* buf)
{
	int ofs = blockDim.x*blockIdx.x   + threadIdx.x;
	float_type maxro_ = calculate_maxplasmadensity_small_device((float_type*)(input+cuda_N_T*ofs), (float_type*)(buf+cuda_N_T*ofs)); 
	output[ofs] = maxro_;
}


#ifdef YUDIN_IVANOV_CORRECTION

__device__ float_type YI_Phi_device(float_type theta, float_type g)
{
 float_type sth = fabs(sin(theta)); 
 float_type sth2 = sth*sth; 

 float_type a = 1.0 + g*g - sth2;
 float_type b = sqrt(a*a + 4*g*g*sth2); 
 float_type c = sqrt(pow(sqrt((a+b)/2.0)+g, 2.0) + pow(sqrt((b-a)/2.0)+sth, 2.0)); 
 
 return ((g*g + sth2 + 0.5)*log(c) - (3.0*sqrt((b-a)/2.0)/2.0)*sth - sqrt((a+b)/2.0)/2.0*g); 
}

#endif



#endif 
