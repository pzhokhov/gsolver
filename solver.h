#ifndef _SOLVER_HEADER
#define _SOLVER_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <malloc.h>

#include <mpi.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include "SIconstants.h"
#include "extmath.h"
#include "dumps.h"
#include "fhatha.h"

#ifdef _ENABLE_CUDA
#include <cuda.h>
#ifndef __float128
	#define  __float128 long long double
#endif
//#include <cutil.h>
#endif

#include <fftw3.h>
#include <fftw3-mpi.h>



#ifdef PPT_IONIZATION
#define IONIZATION_GAS
#endif


//#define NO_PLASMARESPONSE

//#undef TUNNEL_IONIZATION
//#define  MULTIPHOTON_IONIZATION
//#undef KELDYSH_IONIZATION
//#undef  DIRECT_IONIZATION_RATE

//#define _SHOW_EVERY_STEP

#define NONLINEARITY_ON 
//#define _SHOW_EVERY_STEP
//#define NO_AVALANCHE
//#define NO_PLASMARESPONSE
//#define NO_PONDEROMOTIVE_POTENTIAL

//#define THIRD_HARMONICS
//#define _UNWRAP_FREQUENCIES



#ifndef _ODE_EULER	
#ifndef _ODE_HEUN
#ifndef _ODE_RK4
#define _ODE_RK4
#endif
#endif
#endif

#define  MINSTEP_RATIO  ((float_type)1e-10)

#ifndef MAX_TOLERANCE_DEFAULT
 #define  MAX_TOLERANCE_DEFAULT  ((float_type)1e-4)
#endif

#define  IONIZATION_I_LN_TOLERANCE ((float_type)0.003304003148325)
#define  IONIZATION_MIN_I_LN    ((float_type)29.933606208922594)        
#define  IONIZATION_MIN_I       ((float_type)1e13)
#define  IONIZATION_MAX_I_LN    ((float_type)57)
#define  IONIZATION_N           8192



#ifndef INTENSITY_DENOM
 #define INTENSITY_DENOM     1e16
 #define FIELD_DENOM         1e8
 #define INTENSITY_DENOM_LN  ((float_type)36.841361487904734)
#endif


#ifndef IONIZATION_GAS
#define IONRATE_DENOM  ((float_type)1e20)
#define IONRATE_DENOM_LN ((float_type)46.0517018)
#else
#define IONRATE_DENOM    ((float_type)1e10)
#define IONRATE_DENOM_LN ((float_type)23.025850929940457)
#endif

/*
#define INTENSITY_DENOM     ((float_type)1.0)
#define FIELD_DENOM         ((float_type)1.0)
#define INTENSITY_DENOM_LN  ((float_type)0.0)
*/

#define DEFAULT_RESIZE_MAX_TOLERANCE  1e-4
#define DEFAULT_RESIZE_MIN_TOLERANCE  1e-7
#define DEFAULT_RESIZE_FREQUENCY	  10

#define  ABSORBTION_LAYER_WIDTH 0.1
#define  ABSORBTION_LAYER_BETA  10.0



#define MAX_KT2	((float_type)0.25)

extern int PROCESS_N;
extern int PROCESS_RANK;

#ifdef _ENABLE_CUDA
extern MPI_Comm DEVICE_COMM;
#endif 

extern int N_T;			extern float_type  TMIN, TMAX, TSTEP;
extern int N_X;			extern float_type  XMIN, XMAX, XSTEP;
extern int N_Y;			extern float_type  YMIN, YMAX, YSTEP;
extern int N_Z, n_Z;	extern float_type* ZNET;

extern int APERTURE_N; 
extern float_type* APERTURE_Z;
extern float_type* APERTURE_R;

extern dump_manager_class* DUMPMAN;

extern FILE* SC_FID; 

extern double TIME_START;

extern f_complex* WAVENUMBER;
extern f_complex  WAVENUMBER0;
extern float_type* HO_DISPERSION; 
extern float_type *OMEGA;
extern float_type  OMEGA0;
extern float_type OMEGA_MAX, OMEGA_MIN;
extern float_type  GROUP_VELOCITY;
extern f_complex *RAMAN_FUNCTION;
extern f_complex PLASMA_FACTOR;

extern float_type ZSTEP, LIN_ZSTEP;
extern float_type MAX_TOLERANCE;

extern float_type RESIZE_MAX_TOLERANCE;
extern float_type RESIZE_MIN_TOLERANCE;
extern int RESIZE_FREQUENCY;


extern ptrdiff_t MY_NX_FT, MY_NXstart_FT;
extern ptrdiff_t MY_NX, MY_NXstart;
extern ptrdiff_t MY_NY, MY_NYstart;
extern ptrdiff_t MY_SIZE;


extern bool STOP_SIGNAL;

extern unsigned int STEP_N;
extern float_type CURRENT_Z, FINAL_Z;

extern f_complex *FIELD, *BIGBUFFER1, *BIGBUFFER2;

extern f_complex *NL_OUTPUT1, *NL_OUTPUT2, *NL_OUTPUT3, *NL_OUTPUT4, *NL_OUTPUT5, *NL_OUTPUT6;
extern f_complex *NL_INPUT;
extern f_complex *FIELD_REAL_SMALL, *NL_SMALLBUFFER1 , *NL_SMALLBUFFER2, *NL_SMALLBUFFER3, *NL_SMALLBUFFER4; 

extern float_type NONLIN_REFRINDEX;
extern float_type NONLIN_REFRINDEX4;

extern float_type *KERR_PROFILE;
extern float_type *KERR_TH_PROFILE;

extern float_type RAMAN_FRACTION, TAU_RAMAN, OMEGA_RAMAN;
extern float_type AMBIENT_CARRIER_DENSITY;

extern float_type TH_FACTOR;

extern float_type IONIZATION_POTENTIAL;

#ifdef MULTI_LEVEL_IONIZATION
 extern float_type* IONIZATION_POTENTIALS;
 extern float_type* ION_DENSITIES_BUFFER;
 extern float_type* IONIZATION_RATES_BUFFER;
 extern int	    IONIZATION_LEVEL_N;
extern float_type* N2_IONFACTOR;
#endif 


#ifndef MULTI_LEVEL_IONIZATION
extern float_type IONIZATION_RATE_LN[];
extern float_type IONIZATION_I_LN   [];
#else
extern float_type* IONIZATION_RATE_LN;
extern float_type* IONIZATION_I_LN;
#endif

#ifdef MULTIPHOTON_IONIZATION
	extern float_type BETA_MPI_LN;
	extern int    K_MPI;
#endif

#ifdef TUNNEL_IONIZATION
	extern float_type TUNNELING_FIELD;
#endif 

extern float_type AVALANCHE_CROSSSECTION;
extern float_type PONDEROMOTIVE_COEFFICIENT;
extern float_type RECOMBINATION_TAU;
extern float_type COLLISION_TAU;
extern float_type NEUTRAL_DENSITY;
extern f_complex* PLASMA_FUNC;


extern int _ARGC;
extern char** _ARGV;

extern fftwt_plan        FFT_FWPLAN_T,  FFT_BWPLAN_T, FFT_ALLBWPLAN_T;  // local computed fft's
extern fftwt_plan  FFT_FWPLAN_XY, FFT_BWPLAN_XY; // collectively computed fft's

#ifdef QDHT
extern qdht_plan*  HT_PLAN; 
#else
extern fhatha_plan* HT_PLAN;
#endif

#ifdef BLOCH_RESPONSE
extern int  BLOCH_N;
extern float_type BLOCH_TEMP;

extern float_type *BLOCH_ENERGY;
extern float_type *BLOCH_DIPOLE;
extern float_type *BLOCH_G1;
extern float_type *BLOCH_G2; 
#endif


void print_helpmsg();
void process_input(int argc, char** argv);
void print_variables();

int      load_namedint           (FILE* fid, const char* name,                       bool defvalue_present=false, int defvalue=0);
int      load_namednumint        (FILE* fid, const char* name, int num,              bool defvalue_present=false, int defvalue=0);
float_type   load_namedfloat     (FILE* fid, const char* name,                       bool defvalue_present=false, float_type defvalue=0);
float_type   load_namednumfloat  (FILE* fid, const char* name, int num,              bool defvalue_present=false, float_type defvalue=0);
float_type   load_named2numfloat (FILE* fid, const char* name, int num1, int num2,   bool defvalue_present=false, float_type defvalue=0);
void         load_namedstringn   (FILE* fid, const char* name, char* output, int N,  bool defvalue_present=false, const char* defvalue=NULL);
void         load_namednumstringn(FILE* fid, const char* name, char* output, int num,   int N, bool defvalue_present=false, const char* defvalue=NULL);

void load_info(FILE*);
void create_mystartcondition(FILE* fid);
void create_net(float_type Xmin, float_type Xmax, int N, char* nettype, float_type* net);

void load_nonlindata(FILE* fid);

void  dump_mypiece	  (f_complex* piece, char* prefix, int N);
void  dump_mypiece	  (f_complex* piece, char* prefix, float_type D);
void _dump_piece(f_complex* piece, char* filename);
void  dump_centralsectiony_init(char* filename);
void  dump_centralsectiony();

void initialize_variables();
void destroy_variables();
void cuda_load_const();
void cuda_free_const();

void computationcycle();

void propagator();

void calc_ionization_rate();

void calculate_plasmadensity_losses_small(f_complex* field, float_type* pro, int stride=1, f_complex* loss=NULL, f_complex* n2factor=NULL);
void calculate_plasmadensity_2float(f_complex* input, float_type* output, size_t N);
void calculate_maxplasmadensity_2float(f_complex* input, float_type* output, size_t N);

float_type photoionization_rate_ln(float_type I);

void master_primary_zstep();


   #define min(a,b) ((a)<(b)?(a):(b))
   #define max(a,b) ((a)>(b)?(a):(b))

#define ISMASTER (PROCESS_RANK == 0)
#define sign(x)  (((x)>=0)?(1):(-1)) 

    
#ifdef _WIN32 
  
    #define WIN32_LEAN_AND_MEAN

#endif

    

inline int    global_nx(int nx)          {return nx;}
inline int    global_ny(int ny)          {return ny+MY_NYstart;}

inline float_type abs2(float_type re, float_type im) {return re*re + im*im;}
inline float_type abs2(float_type* z)            {return (*z)*(*z) + (*(z+1))*(*(z+1));}
inline float_type abs2(f_complex z)     {return real(z)*real(z)+imag(z)*imag(z);}

inline void  fftwt_Nnormalize(int pN, float_type* M, int N=N_T)         {for(long i=0; i<2*pN*N; i++) {M[i]/=N;}}
inline void  fftwt_Nnormalize(int pN, f_complex* M, int N=N_T) {for(long i=0; i<  pN*N; i++) {M[i]/=(float_type)N;}}

inline float_type getmaxI(f_complex* A, int N) {float_type maxI = abs2(A[0]); for(int i=1; i<N; i++) maxI=max(maxI, abs2(A[i])); return maxI;}
inline float_type getminI(f_complex* A, int N) {float_type minI = abs2(A[0]); for(int i=1; i<N; i++) minI=min(minI, abs2(A[i])); return minI;}


inline void fill_zeros(float_type* M,         size_t N) {for (size_t i=0; i < N; i++) M[i] = 0.0;}
inline void fill_zeros(f_complex* M, size_t N) {for (size_t i=0; i < N; i++) M[i] = 0.0;}

inline float_type hfgaussfilter(float_type omega) {float_type domega = (omega/OMEGA_MAX-1)/ABSORBTION_LAYER_WIDTH; return (domega<0)?(float_type)1.0:exp(-ABSORBTION_LAYER_BETA*domega*domega);} 
inline float_type lfgaussfilter(float_type omega) {float_type domega = (omega*(1-ABSORBTION_LAYER_WIDTH)/OMEGA_MIN-1); return (domega>0)?(float_type)1.0:exp(-ABSORBTION_LAYER_BETA*domega*domega);} 

//inline float_type hfgaussfilter(float_type omega) {float_type domega = (2*omega/OMEGA_MAX-1)/0.99; return (float_type)exp(-pow(domega, 100));} 
//inline float_type lfgaussfilter(float_type omega) {return 1;} 


inline float_type frand_common() {float_type f=0; if (ISMASTER) f=frand(); MPI_Bcast(&f, 1, MPI_FLOAT_TYPE, 0, MPI_COMM_WORLD); return f;}

#ifdef QDHT
  #define hankel_runmany_cuda      qdht_runmany_cuda
  #define hankel_runmany_MPI       qdht_runmany_MPI
  #define hankel_runmany_MPI_cuda  qdht_runmany_MPI_cuda
#else 
  #define hankel_runmany_cuda      fhatha_runmany_cuda
  #define hankel_runmany_MPI       fhatha_runmany_MPI
  #define hankel_runmany_MPI_cuda  fhatha_runmany_MPI_cuda
#endif

#ifdef _ENABLE_CUDA 
  extern bool cuda_do_extra_memtransfer; 
#endif

inline void ht_run(f_complex* data, f_complex* buf) 
{
#ifdef _ENABLE_CUDA
	 if (PROCESS_N == 1) hankel_runmany_cuda(HT_PLAN, data, N_T, N_T, 1, cuda_do_extra_memtransfer, (float_type*)buf);
	 else hankel_runmany_MPI_cuda(HT_PLAN, data, N_T, N_T, 1, buf, MPI_COMM_WORLD);
#else
	 if (PROCESS_N == 1) HT_PLAN->run_many(data, N_T, N_T, 1); 
	 else hankel_runmany_MPI(HT_PLAN, data, N_T, N_T, 1, buf, MPI_COMM_WORLD);
#endif 
}


void load_spectrum_fromfile(f_complex* output,float_type* omega, int Nt,  const char* filename, const char* filetype);

#ifdef _UNIAXIAL_FINITE_DIFFERENCE
 void calculate_Lstep_UAFD();
#endif

#ifndef TRANSVERSE_DIMENSIONS
  #error Number of transverse dimensions should be specified in TRANSVERSE_DIMENSIONS macro
#endif


#define printf_fl(...) {printf(__VA_ARGS); fflush(stdout);}

#endif
