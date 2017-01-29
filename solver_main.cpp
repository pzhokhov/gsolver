// solver_main.cpp : Defines the entry point for the console application.
// The main, parallelized program

#include "solver.h"
#include "ionization.h"

int PROCESS_N;
int PROCESS_RANK;
MPI_Comm DEVICE_COMM;

int N_T;		float_type  TMIN, TMAX, TSTEP;
int N_X;		float_type  XMIN, XMAX, XSTEP;
int N_Y;		float_type  YMIN, YMAX, YSTEP;
int N_Z, n_Z=0;	float_type* ZNET;

int APERTURE_N=0; float_type* APERTURE_Z=NULL; float_type* APERTURE_R=NULL; 

FILE*  CSDUMP_FID = NULL, *CSDUMP_PLASMA_FID = NULL;

FILE*  SC_FID = NULL;

dump_manager_class* DUMPMAN = NULL;

double TIME_START = 0;

float_type LAMBDA_V   = 0;

float_type ZSTEP = 0, LIN_ZSTEP = 0;
float_type MAX_TOLERANCE = MAX_TOLERANCE_DEFAULT;


float_type RESIZE_MAX_TOLERANCE = DEFAULT_RESIZE_MAX_TOLERANCE;
float_type RESIZE_MIN_TOLERANCE = DEFAULT_RESIZE_MIN_TOLERANCE;
int RESIZE_FREQUENCY = DEFAULT_RESIZE_FREQUENCY;

ptrdiff_t MY_NX_FT, MY_NXstart_FT;
ptrdiff_t MY_NX, MY_NXstart; 
ptrdiff_t MY_NY, MY_NYstart;
ptrdiff_t MY_SIZE;


bool STOP_SIGNAL = false;

unsigned int STEP_N;
float_type CURRENT_Z = 0, FINAL_Z = 0;

f_complex *FIELD, *BIGBUFFER1, *BIGBUFFER2;

f_complex *NL_OUTPUT1, *NL_OUTPUT2, *NL_OUTPUT3, *NL_OUTPUT4, *NL_OUTPUT5, *NL_OUTPUT6;
f_complex *NL_INPUT;
f_complex *FIELD_REAL_SMALL, *NL_SMALLBUFFER1 , *NL_SMALLBUFFER2, *NL_SMALLBUFFER3, *NL_SMALLBUFFER4;



float_type *OMEGA;
float_type  OMEGA0;
float_type OMEGA_MIN = 0, OMEGA_MAX = 1e20;
f_complex *WAVENUMBER;
f_complex  WAVENUMBER0;
float_type GROUP_VELOCITY;
float_type *HO_DISPERSION;

f_complex *RAMAN_FUNCTION;
f_complex PLASMA_FACTOR;

float_type AMBIENT_CARRIER_DENSITY = 0;


float_type* KERR_PROFILE;

float_type* KERR_TH_PROFILE;


float_type NONLIN_REFRINDEX = 0;
float_type NONLIN_REFRINDEX4 = 0;

float_type RAMAN_FRACTION = 0, TAU_RAMAN = 0, OMEGA_RAMAN = 0;

float_type TH_FACTOR = 0;

#ifndef MULTI_LEVEL_IONIZATION
 float_type IONIZATION_RATE_LN [IONIZATION_N];
 float_type IONIZATION_I_LN    [IONIZATION_N];
#else
 float_type* IONIZATION_RATE_LN ;
 float_type* IONIZATION_I_LN    ;
#endif 
float_type IONIZATION_POTENTIAL;

#ifdef MULTI_LEVEL_IONIZATION
 float_type* IONIZATION_POTENTIALS;
 float_type* ION_DENSITIES_BUFFER;
 float_type* IONIZATION_RATES_BUFFER;
 int	    IONIZATION_LEVEL_N;
 float_type* N2_IONFACTOR;
#endif 


#ifdef MULTIPHOTON_IONIZATION
	float_type BETA_MPI_LN;
	int    K_MPI;
#endif

#ifdef TUNNEL_IONIZATION
	float_type TUNNELING_FIELD;
#endif 


float_type AVALANCHE_CROSSSECTION;
float_type PONDEROMOTIVE_COEFFICIENT;
float_type NEUTRAL_DENSITY;
float_type RECOMBINATION_TAU;
float_type COLLISION_TAU;

f_complex* PLASMA_FUNC;

int _ARGC;
char** _ARGV;

fftwt_plan        FFT_FWPLAN_T,  FFT_BWPLAN_T, FFT_ALLBWPLAN_T;  
fftwt_plan  FFT_FWPLAN_XY, FFT_BWPLAN_XY; 

#if TRANSVERSE_DIMENSIONS == 1
 #ifdef QDHT
  qdht_plan* HT_PLAN;
 #else
  fhatha_plan* HT_PLAN;
 #endif
#endif 


#ifdef BLOCH_RESPONSE
int BLOCH_N;
float_type BLOCH_TEMP;
float_type *BLOCH_RHO_D; 
f_complex  *BLOCH_RHO_O; 
float_type *BLOCH_ENERGY;
float_type *BLOCH_DIPOLE;
float_type *BLOCH_G1;
float_type *BLOCH_G2; 

#endif

void calc_ho_dispersion(); 
void calc_raman_function();

void resize_if_needed();

inline void print_status(double timetot) { 	  if (ISMASTER) {printf("\nStep number %d, total time %5.2fs, time per step %3.2fs, ETA %5.2fs", STEP_N, timetot, timetot/STEP_N, timetot*((double)N_Z/(double)n_Z - 1.0));fflush(stdout);	  }}

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &PROCESS_N);
	MPI_Comm_rank(MPI_COMM_WORLD, &PROCESS_RANK);
	
    _ARGC = argc;
    _ARGV = argv;
	
	fftwt_mpi_init();
                                  	
	if (ISMASTER) 
	{
#ifdef _ENABLE_CUDA
	 printf("*** UPPE solver (CUDA GPU version)*** \n");
#else 
	 printf("*** UPPE solver (CPU version, on %d CPUs)*** \n", PROCESS_N);
#endif

#ifdef  _EXTMATH_SINGLE
	 printf("Float-type precision: single\n");
#else
	 printf("Float-type precision: double\n");
#endif
	 time_t starttime; time(&starttime);
	 printf("\n %s", ctime(&starttime));
#ifdef _FULL_DUMPS		
	 srand((unsigned)time(NULL));
	 DUMP_ID = rand();
	 printf("\n Dump id is %d",DUMP_ID);
	 fflush(stdout);
#endif
	}
#ifdef _SINGLE_PROCESS
	if (PROCESS_N > 1) throw "The source code was compiled with _SINGLE_PROCESS key, it does not support parallel execution";
#endif
	srand48((unsigned)(time(NULL)+PROCESS_RANK));

	MPI_Barrier(MPI_COMM_WORLD);
	try
	{
	 STEP_N = 0;
	 process_input(argc, argv);
	 TIME_START = MPI_Wtime();
	 if (ISMASTER) printf("\n ZSTEP = %e", ZSTEP);

#ifndef _SILENCE
	 if (ISMASTER) printf("\nStarting computation cycle");
#endif
	 double status_time = TIME_START; 
         while (!STOP_SIGNAL)
	 {

	  double timetot = MPI_Wtime()-TIME_START; if (timetot > status_time + 30) {print_status(timetot); status_time = timetot;}
	  if (CURRENT_Z >= ZNET[n_Z])
	  {
		  DUMPMAN -> dump();
		  print_status(timetot); status_time = timetot; 

    	  if (n_Z == N_Z-1) {STOP_SIGNAL = true; break;}
	      n_Z++;
	  }			
			
#ifndef _SILENCE
	  fflush(stdout);
#endif
	  computationcycle();

	  CURRENT_Z += ZSTEP;
#ifdef ONECYCLE
	  STOP_SIGNAL = true;
#endif
				
	 }
#ifdef DUMP_FINAL_FIELD
	 fftw(FFT_BWPLAN_T, N_X*MY_NY, (fftwt_complex*)FIELD, 1,N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(N_X*MY_NY, BIGBUFFER1);
	 dump_mypiece(BIGBUFFER1, DUMP_PREFIX, n_Z);
#endif
	 MPI_Barrier(MPI_COMM_WORLD);
	 fflush(stdout);
	 destroy_variables();	
	}
	catch (char* errmsg)
	{
 	 printf("\n!!!Error in process %d occured!!! , %s", PROCESS_RANK, errmsg);
         fflush(stdout);
	 MPI_Abort(MPI_COMM_WORLD,2);
	 return 1;
	}
	catch (const char* errmsg)
	{
	 printf("\n!!!Error in process %d occured!!! , %s", PROCESS_RANK, errmsg);
	 fflush(stdout);
	 MPI_Abort(MPI_COMM_WORLD,2);
	 return 1;
	}

	MPI_Finalize();
	return 0;
}

void computationcycle()
{
 propagator(); 
 #ifndef NO_RESIZE
 if (STEP_N % RESIZE_FREQUENCY == 0)
 {
	 resize_if_needed();
 }
 #endif
 STEP_N++;
}

void initialize_variables()
{
#if TRANSVERSE_DIMENSIONS == 0
	MY_NX = 1;
	MY_NY = 1;
	MY_NXstart = 0;
	MY_NYstart = 0;
	MY_NX_FT = 1;
	MY_SIZE = N_T; 
	
#elif TRANSVERSE_DIMENSIONS == 1
#ifdef QDHT
	HT_PLAN =  new qdht_plan(N_X);
	if (ISMASTER) printf("\n Hankel transform - quasi-discrete. ");

#else
	float_type fresnel_number = load_namedfloat(SC_FID, "FRESNEL_NUMBER",true,-1);
	if (ISMASTER && fresnel_number==-1) printf("\n Warning!!! Fresnel number is not specified. Default might not be the best one!");

	float_type alpha_mult = load_namedfloat(SC_FID, "R_NET_ALPHA_FACTOR", true, 1);
	HT_PLAN =  new fhatha_plan(N_X, alpha_mult);

	if (ISMASTER) printf("\n Hankel transform parameters:  Fresnel number=%.4g, r0=%3e m, rmax=%3e m", HT_PLAN->getNf(), XMAX*HT_PLAN->x_n(0), XMAX*HT_PLAN->x_n(N_X-1));
#endif


	MY_NX = N_X/PROCESS_N;
	MY_NXstart = MY_NX*PROCESS_RANK;
	MY_SIZE = MY_NX*N_T;
	MY_NY = 1;
	MY_NYstart = 0; 
	MY_NX_FT = MY_NX;
	MY_NXstart_FT = MY_NXstart;
	
#elif TRANSVERSE_DIMENSIONS == 2
	ptrdiff_t xyplan_sizes[2] = {N_Y, N_X}; 	
	MY_SIZE = fftwt_mpi_local_size_many_transposed(2, xyplan_sizes, N_T, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &MY_NY, &MY_NYstart, &MY_NX_FT, &MY_NXstart_FT); 

	MY_NXstart = 0;
	MY_NX = N_X;
#endif


	unsigned long datapiece_size       = sizeof(f_complex)*MY_SIZE;
	unsigned long datapiece_small_size = sizeof(f_complex)*N_T;

    BIGBUFFER1    = (f_complex*)malloc_ch(datapiece_size);  
#ifndef _UNIAXIAL_FINITE_DIFFERENCE
   	FIELD         = (f_complex*)malloc_ch(datapiece_size);  
	BIGBUFFER2    = (f_complex*)malloc_ch(datapiece_size); 
#else
   	FIELD         = (f_complex*)malloc_ch(N_T*(MY_NX+1)*sizeof(f_complex));  
	BIGBUFFER2    = (f_complex*)malloc_ch(2*datapiece_size);
#endif

#if TRANSVERSE_DIMENSIONS == 2
	FFT_FWPLAN_XY = fftwt_mpi_plan_many_dft(2, xyplan_sizes, N_T, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, (fftwt_complex*)BIGBUFFER1,  (fftwt_complex*)BIGBUFFER1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MPI_TRANSPOSED_OUT | FFTW_FLAG);
	FFT_BWPLAN_XY = fftwt_mpi_plan_many_dft(2, xyplan_sizes, N_T, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, (fftwt_complex*)FIELD,  (fftwt_complex*)FIELD, MPI_COMM_WORLD, FFTW_FORWARD , FFTW_MPI_TRANSPOSED_IN  | FFTW_FLAG);
#endif

    FIELD_REAL_SMALL		= (f_complex*)malloc_ch(datapiece_small_size);
	NL_SMALLBUFFER1			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_SMALLBUFFER2			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_SMALLBUFFER3			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_SMALLBUFFER3			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_SMALLBUFFER4			= (f_complex*)malloc_ch(datapiece_small_size);

	RAMAN_FUNCTION			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_OUTPUT1			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_OUTPUT2			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_OUTPUT3			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_OUTPUT4			= (f_complex*)malloc_ch(datapiece_small_size);
	NL_INPUT			= (f_complex*)malloc_ch(datapiece_small_size);

	FFT_FWPLAN_T    = fftwt_plan_dft_1d(N_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)NL_SMALLBUFFER2, FFTW_FORWARD,  FFTW_FLAG);
	FFT_BWPLAN_T    = fftwt_plan_dft_1d(N_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)NL_SMALLBUFFER2, FFTW_BACKWARD, FFTW_FLAG);
	FFT_ALLBWPLAN_T = fftwt_plan_many_dft(1, &N_T, MY_NX*MY_NY, (fftwt_complex*)FIELD, NULL, 1, N_T, (fftwt_complex*)BIGBUFFER1,NULL,  1, N_T, FFTW_BACKWARD, FFTW_FLAG);
	
	KERR_PROFILE			= (float_type*)malloc_ch(datapiece_small_size);
	KERR_TH_PROFILE		        = (float_type*)malloc_ch(datapiece_small_size);

	CURRENT_Z = ZNET[0];
	
	AVALANCHE_CROSSSECTION = 0;
	PONDEROMOTIVE_COEFFICIENT = 0;

	f_complex j = f_complex(0,1);
	float_type m = ELECTRON_MASS*load_namedfloat(SC_FID, "REDUCED_MASS", true, 1); 
	PLASMA_FACTOR = (ELECTRON_CHARGE*ELECTRON_CHARGE*VACUUM_PERMEABILITY/(WAVENUMBER0*WAVENUMBER0*m*((f_complex)1.0-j/OMEGA0/COLLISION_TAU)));

	float_type IBS_crossection = ELECTRON_CHARGE/m*ELECTRON_CHARGE*VACUUM_PERMEABILITY/real(WAVENUMBER0);
#ifndef NO_AVALANCHE
	
	AVALANCHE_CROSSSECTION =    IBS_crossection*OMEGA0*COLLISION_TAU/(1+OMEGA0*COLLISION_TAU*OMEGA0*COLLISION_TAU);

#endif
#ifndef NO_PONDEROMOTIVE_POTENTIAL
    PONDEROMOTIVE_COEFFICIENT = IBS_crossection/2.0/real(WAVENUMBER0)/GROUP_VELOCITY/(1.0+1.0/(OMEGA0*COLLISION_TAU));
#endif

	PLASMA_FUNC = (f_complex*)malloc(sizeof(f_complex)*N_T);
	for (int nw=0; nw<N_T; nw++)
	{
		float_type alpha =  COLLISION_TAU*OMEGA[nw];
		float_type alpha0 = COLLISION_TAU*OMEGA0;
		float_type k = real(WAVENUMBER[nw]);

		//alpha = COLLISION_TAU*3e14; 
		//k = 3e14/LIGHT_VELOCITY;
		PLASMA_FUNC[nw] = -(alpha/((float_type)1.0+j*alpha)/(float_type)2.0/k*ELECTRON_CHARGE*ELECTRON_CHARGE/m*VACUUM_PERMEABILITY);
		
		//PLASMA_FUNC[nw] = -alpha0/((float_type)1.0+j*alpha)*IBS_crossection/(float_type)2.0;
#ifdef PLASMA_FULL_DISPERSION
		PLASMA_FUNC[nw] *= (float_type)2.0*j/k;
#endif
	}

	calc_ionization_rate();
	calc_raman_function();         

	IONIZATION_POTENTIAL *=IONRATE_DENOM/INTENSITY_DENOM;
#ifdef MULTI_LEVEL_IONIZATION
	for (size_t i=0; i<IONIZATION_LEVEL_N; i++) IONIZATION_POTENTIALS[i]*=IONRATE_DENOM/INTENSITY_DENOM;
#endif

	NONLIN_REFRINDEX*=INTENSITY_DENOM;
	NONLIN_REFRINDEX4*=INTENSITY_DENOM*INTENSITY_DENOM/1e16/1e16;
	RECOMBINATION_TAU *= IONRATE_DENOM;
	PONDEROMOTIVE_COEFFICIENT *= IONRATE_DENOM;

#ifdef _ENABLE_CUDA
	cuda_load_const();
#endif 
}


void calc_raman_function()
{
	float_type M = TSTEP*(1+OMEGA_RAMAN*OMEGA_RAMAN*TAU_RAMAN*TAU_RAMAN)/TAU_RAMAN/(OMEGA_RAMAN*TAU_RAMAN);
	for (int nt=0; nt<N_T; nt++)
	{
		float_type t = TSTEP*nt;
		NL_SMALLBUFFER1[nt] = RAMAN_FRACTION*M*exp(-t/TAU_RAMAN)*sin(OMEGA_RAMAN*t);

	}
	fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)RAMAN_FUNCTION);
	//fftwt_Nnormalize(1, RAMAN_FUNCTION); 
	for (int nw=0; nw<N_T; nw++) RAMAN_FUNCTION[nw] += (1-RAMAN_FRACTION);
}

void destroy_variables()
{
	fclose(SC_FID);
#ifdef _ENABLE_CUDA
	cuda_free_const();
#endif

        free(FIELD);
  	free(BIGBUFFER1);
        free(BIGBUFFER2);

	free(OMEGA); free(WAVENUMBER); free(HO_DISPERSION);

	free(FIELD_REAL_SMALL); free(NL_SMALLBUFFER1);free(NL_SMALLBUFFER2);free(NL_SMALLBUFFER3); free(NL_SMALLBUFFER4);
	free(KERR_PROFILE);     free(KERR_TH_PROFILE);

	free(RAMAN_FUNCTION); free(PLASMA_FUNC); free(APERTURE_R); free(APERTURE_Z);

#ifdef MULTI_LEVEL_IONIZATION
	free(IONIZATION_POTENTIALS); free(ION_DENSITIES_BUFFER); free(IONIZATION_RATES_BUFFER); free(N2_IONFACTOR);
#endif

	fftwt_destroy_plan(FFT_FWPLAN_T); 	fftwt_destroy_plan(FFT_BWPLAN_T);
	fftwt_destroy_plan(FFT_FWPLAN_XY);       fftwt_destroy_plan(FFT_BWPLAN_XY);
	delete DUMPMAN;
}


void resize_if_needed()
{
/*	float_type maxIt = 0,  maxIw = 0,  maxpIt = 0,  maxpIw = 0,  maxpIt2 = 0,  maxpIw2 = 0; 
	float_type maxIt_ = 0, maxIw_ = 0, maxpIt_ = 0, maxpIw_ = 0, maxpIt2_ = 0, maxpIw2_ = 0; 

	float_type maxpIx  = 0,  maxpIkx = 0,  maxpIx2  = 0, maxpIkx2 = 0; 
	float_type maxpIx_ = 0,  maxpIkx_ = 0, maxpIx2_ = 0, maxpIkx2_ = 0; 

	float_type maxpIy  = 0,  maxpIky = 0,  maxpIy2  = 0, maxpIky2 = 0; 
	float_type maxpIy_ = 0,  maxpIky_ = 0, maxpIy2_ = 0, maxpIky2_ = 0; 

	bool resized = false; 

	fftw(FFT_BWPLAN_T, MY_NX*MY_NY, (fftwt_complex*)FIELD, 1,N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(MY_NX*MY_NY, BIGBUFFER1);

	if (maxpIt > RESIZE_MAX_TOLERANCE*maxIt && maxpIw2 < RESIZE_MIN_TOLERANCE*maxIw) 
	{
		expandtwice(BIGBUFFER1, 1, N_T, N_T, MY_NX*MY_NY, 2);  // expand field in t domain
		fftw(FFT_FWPLAN_T, MY_NX*MY_NY, (fftwt_complex*)BIGBUFFER1, 1,N_T, (fftwt_complex*)FIELD, 1,N_T); 

		fftw(FFT_BWPLAN_T, MY_NX*MY_NY, (fftwt_complex*)BIGBUFFER2, 1,N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(MY_NX*MY_NY, BIGBUFFER1);
		expandtwice(BIGBUFFER1, 1, N_T, N_T, MY_NX*MY_NY, 2);  // expand x-y Fourier (or Hankel) transform in temporal domain
		fftw(FFT_FWPLAN_T, MY_NX*MY_NY, (fftwt_complex*)BIGBUFFER1, 1,N_T, (fftwt_complex*)BIGBUFFER2, 1,N_T); 
		T_MIN = T_MIN*2;
		T_MAX = T_MAX*2;
		recalculate_profiles();
	}

	if (maxpIw > RESIZE_MAX_TOLERANCE*maxIw && maxpIt2 < RESIZE_MIN_TOLERANCE*maxIt) 
	{
		expandtwice(FIELD,      1, N_T, N_T, MY_NX*MY_NY, 1);  // expand field in FT_t domain
		expandtwice(BIGBUFFER2, 1, N_T, N_T, MY_NX*MY_NY, 1);  // expand x-y Fourier (or Hankel) transform in FT_t domain
		T_MIN = T_MIN / 2;
		T_MAX = T_MAX / 2;
		recalculate_profiles();
	}
	
	for (int i=0; i<N_T*MY_NX*MY_NY; i++) {maxIt = max(maxIt, abs2(FIELD[i])); maxIw = max(maxIw, abs2(BIGBUFFER2[i]));}
	maxIt_ = maxIt; maxIw_ = maxIw;
	
	MPI_Allreduce(&maxIt_,   &maxIt,   1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  
	MPI_Allreduce(&maxIw_,   &maxIw,   1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  

#ifndef _UNIAXIAL
	//
	//x-resizing


	
	for (int nt=0; nt<N_T; nt++)
	{
		for (int ny=0; ny<MY_NY; ny++) 
		{
			int ofs0 = nt+N_T*N_X*ny; 
			maxpIx   = max(maxpIx,  abs2(FIELD[ofs0]));				   maxpIx   = max(maxpIx,  abs2(FIELD[ofs0+N_T*(N_X-1)]));
			maxpIx2  = max(maxpIx2, abs2(FIELD[ofs0 + N_T*3*N_X/4]));  maxpIx2  = max(maxpIx2, abs2(FIELD[ofs0 + N_T*N_X/4])); 
		}
		for (int nx=0; nx<MY_NX_FT; nx++)
		{
			int ofs0 = nt+N_T*N_Y*nx; 
			maxpIky   = max(maxpIky,  abs2(BIGBUFFER2[ofs0]));				  maxpIky   = max(maxpIky,  abs2(BIGBUFFER2[ofs0+N_T*(N_Y-1)]));
			maxpIky2  = max(maxpIky2, abs2(BIGBUFFER2[ofs0 + N_T*3*N_Y/4]));  maxpIky2  = max(maxpIky2, abs2(BIGBUFFER2[ofs0+N_T*N_Y/4])); 
		}
	}

	transpose_mpi(FFT_FWPLAN_XY->p_transpose_inv, sizeof(f_complex)*N_T/sizeof(TRANSPOSE_EL_TYPE), (TRANSPOSE_EL_TYPE*)BIGBUFFER2, (TRANSPOSE_EL_TYPE*)BIGBUFFER1); 
	
	for (int nt=0; nt<N_T; nt++)
	{
		for (int ny=0; ny<MY_NY; ny++) 
		{
			int ofs0 = nt+N_T*N_X*ny; 
			maxpIkx   = max(maxpIkx,  abs2(BIGBUFFER2[ofs0]));				   maxpIkx   = max(maxpIkx,  abs2(BIGBUFFER2[ofs0+N_T*(N_X-1)]));
			maxpIkx2  = max(maxpIkx2, abs2(BIGBUFFER2[ofs0 + N_T*3*N_X/4]));   maxpIkx2  = max(maxpIkx2, abs2(BIGBUFFER2[ofs0 + N_T*N_X/4])); 
		}
	}
	

	maxpIx_ = maxpIx; maxpIx2_ = maxpIx2;  maxpIkx_ = maxpIkx; maxpIkx2_ = maxpIkx2;
	
	MPI_Allreduce(&maxpIx_,   &maxpIx,   1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  
	MPI_Allreduce(&maxpIx2_,  &maxpIx2,  1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  
	MPI_Allreduce(&maxpIkx_,  &maxpIkx,  1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  
	MPI_Allreduce(&maxpIkx2_, &maxpIkx2, 1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  

	if (maxpIx > RESIZE_MAX_TOLERANCE*maxIt && maxpIkx2 < RESIZE_MIN_TOLERANCE*maxIw) 
	{
		for (int ny=0; ny<MY_NY; ny++) expandtwice(FIELD+N_T*N_X*ny, N_T, 1, N_X, N_T, 2);  // expand field in x domain
		float_type xmin = XMIN, xmax = XMAX; 
		XMIN = (xmin + xmax)/2 - (xmax-xmin);
		XMAX = (xmin + xmax)/2 + (xmax-xmin);

		if (ISMASTER) printf("\n X grid is twice expanded.");
		resized = true;
	}

	if (maxpIkx > RESIZE_MAX_TOLERANCE*maxIw && maxpIx2 < RESIZE_MIN_TOLERANCE*maxIt) 
	{
		for (int ny=0; ny<MY_NY; ny++) contracttwice(FIELD+N_T*N_X*ny, N_T, 1, N_X, N_T, 2);  // contract field in x domain
		float_type xmin = XMIN, xmax = XMAX; 
		XMIN = (xmin + xmax)/2 - (xmax-xmin)/4;
		XMAX = (xmin + xmax)/2 + (xmax-xmin)/4;

		resized = true; 
		if (ISMASTER) printf("\n X grid is twice contracted.");
	}
	

	//y-resizing

	transpose_mpi(FFT_FWPLAN_XY->p_transpose, sizeof(f_complex)*N_T/sizeof(TRANSPOSE_EL_TYPE), (TRANSPOSE_EL_TYPE*)FIELD, (TRANSPOSE_EL_TYPE*)BIGBUFFER1); 
	
	for (int nt=0; nt<N_T; nt++)
	{
		for (int nx=0; nx<MY_NX_FT; nx++) 
		{
			int ofs0 = nt+N_T*N_Y*nx; 
			maxpIy   = max(maxpIy,  abs2(FIELD[ofs0]));				   maxpIy   = max(maxpIy,  abs2(FIELD[ofs0 + N_T*(N_Y-1)]));
			maxpIy2  = max(maxpIy2, abs2(FIELD[ofs0 + N_T*3*N_Y/4]));  maxpIy2  = max(maxpIy2, abs2(FIELD[ofs0 + N_T* N_Y/4])); 
		}
	}

	maxpIy_ = maxpIy; maxpIy2_ = maxpIy2;  maxpIky_ = maxpIky; maxpIky2_ = maxpIky2;
	MPI_Allreduce(&maxpIy_,   &maxpIy,   1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  
	MPI_Allreduce(&maxpIy2_,  &maxpIy2,  1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  
	MPI_Allreduce(&maxpIky_,  &maxpIky,  1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  
	MPI_Allreduce(&maxpIky2_, &maxpIky2, 1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);

	if (maxpIy > RESIZE_MAX_TOLERANCE*maxIt && maxpIky2 < RESIZE_MIN_TOLERANCE*maxIw) 
	{
		for (int nx=0; nx<MY_NX_FT; nx++) expandtwice(FIELD+N_T*N_Y*nx, N_T, 1, N_Y, N_T, 2);  // expand field in y domain
		float_type ymin = YMIN, ymax = YMAX; 
		YMIN = (ymin + ymax)/2 - (ymax-ymin);
		YMAX = (ymin + ymax)/2 + (ymax-ymin);


		resized = true; 
		if (ISMASTER) printf("\n Y grid is twice expanded.");
	}

	
	if (maxpIky > RESIZE_MAX_TOLERANCE*maxIw && maxpIy2 < RESIZE_MIN_TOLERANCE*maxIt) 
	{
		for (int nx=0; nx<MY_NX_FT; nx++) contracttwice(FIELD+N_T*N_Y*nx, N_T, 1, N_Y, N_T, 2);  // contract field in y domain
		float_type ymin = YMIN, ymax = YMAX; 
		YMIN = (ymin + ymax)/2 - (ymax-ymin)/4;
		YMAX = (ymin + ymax)/2 + (ymax-ymin)/4;

		resized = true; 
		if (ISMASTER) printf("\n Y grid is twice contracted.");
	}
	
	transpose_mpi(FFT_FWPLAN_XY->p_transpose_inv, sizeof(f_complex)*N_T/sizeof(TRANSPOSE_EL_TYPE), (TRANSPOSE_EL_TYPE*)FIELD, (TRANSPOSE_EL_TYPE*)BIGBUFFER1); 
	
	if (resized) {memcpy(BIGBUFFER2, FIELD, sizeof(f_complex)*N_T*MY_NX*MY_NY); fftwnd_mpi(FFT_FWPLAN_XY, N_T, (fftwt_complex*)BIGBUFFER2, (fftwt_complex*)BIGBUFFER1, FFTW_TRANSPOSED_ORDER);}
	else transpose_mpi(FFT_FWPLAN_XY->p_transpose, sizeof(f_complex)*N_T/sizeof(TRANSPOSE_EL_TYPE), (TRANSPOSE_EL_TYPE*)BIGBUFFER2, (TRANSPOSE_EL_TYPE*)BIGBUFFER1);  

#endif */
}


