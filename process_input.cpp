#include "solver.h"
#include "pulseshapes.h"

#include <gsl/gsl_interp.h>

void calculate_omega(int Nt, float_type tmin, float_type tmax, float_type omega0, float_type* omega);
void calculate_wavenumber(char* filename, int N, float_type* omega, f_complex* wavenum);

float_type calculate_groupvelocity(int Nt, float_type* omega, f_complex* wavenum, float_type omega0);

void calculate_HOdispersion();

void load_refrindex_sellmeier_omega (FILE* fid, int N, float_type* cfreq, f_complex* refrindex);
void load_refrindex_sellmeier_lambda(FILE* fid, int N, float_type* cfreq, f_complex* refrindex);
void load_refrindex_raw             (FILE* fid, int N, float_type* cfreq, f_complex* refrindex);

void load_dispersion_row(FILE* fid, int N, float_type* omega, f_complex* wavenum);

float_type integrate_gaussian(float_type xmin, float_type xmax, float_type width);

void load_kerr_th_profiles(FILE* scfid);

void init_zstep_kerr();

enum 		{FILETYPE_NOFILE           = 0,  \
        	 FILETYPE_SELLMEIER_LAMBDA = 1,  \
		 FILETYPE_SELLMEIER_OMEGA  = 2,  \
		 FILETYPE_RAW              = 3,  \
		 FILETYPE_DISPROW   	   = 10, \
		 FILETYPE_KERR_PROFILE     = 11, \
		 FILETYPE_KERR_TH_PROFILE  = 12};

void process_input(int argc, char** argv)
{
	if (argc < 2 && PROCESS_RANK == 0) throw "\n!!! Error: Too few arguements  !!!";	

	#ifndef _SILENCE
	if (ISMASTER) printf("\n Loading starting info from %s", argv[1]);
	#endif
	   
	SC_FID = fopen(argv[1], "r");
	if (SC_FID == NULL) throw "Error opening file with starting condition.";
	load_info(SC_FID);
	load_nonlindata(SC_FID);

	#ifndef _SILENCE
	if (ISMASTER) printf("\n Initializing variables...");
	#endif
	initialize_variables();
	load_kerr_th_profiles(SC_FID); 

    if (ISMASTER) print_variables();	

	DUMPMAN = new dump_manager_class(SC_FID);
	
	create_mystartcondition(SC_FID);	

if (INTENSITY_DENOM != 1) for (size_t nx=0; nx<MY_NX; nx++) for (size_t ny=0; ny<MY_NY; ny++) for (size_t nw=0; nw<N_T; nw++) FIELD[nw+N_T*(nx+MY_NX*ny)]/=FIELD_DENOM;
    
	
	// create_mystartcondition() outputs Fourier-transform of the field to FIELD array. 
	// To calculate initial zstep, one needs to perform an inverse transform to BIGBUFFER1
	//
	
	fftwt_execute(FFT_ALLBWPLAN_T); fftwt_Nnormalize(MY_NX*MY_NY, BIGBUFFER1);
	init_zstep_kerr();

	memcpy(BIGBUFFER2, FIELD, sizeof(f_complex)*N_T*MY_NX*MY_NY);

#ifndef _UNIAXIAL_FINITE_DIFFERENCE

#if   TRANSVERSE_DIMENSIONS == 1
	ht_run(BIGBUFFER2, BIGBUFFER1);
#elif TRANSVERSE_DIMENSIONS == 2
	fftwt_mpi_execute_dft(FFT_FWPLAN_XY, (fftwt_complex*)BIGBUFFER2, (fftwt_complex*)BIGBUFFER2);
#endif

#else
	/*for (int x=0; x==0;)
 	{
	 x =0;
	}*/
	MPI_Barrier(MPI_COMM_WORLD);
	float_type h = XMAX/N_X;
	LIN_ZSTEP = MAX_TOLERANCE*h*h*real(WAVENUMBER[0]);
	for (int nw=0; nw<N_T; nw++)
	{
	  MPI_Status mpistat; 
	  if (PROCESS_RANK < PROCESS_N-1) MPI_Recv(FIELD+nw+MY_NX*N_T, 2, MPI_FLOAT_TYPE, PROCESS_RANK+1, nw, MPI_COMM_WORLD, &mpistat);
	  else FIELD[nw+MY_NX*N_T]=0;
	  if (PROCESS_RANK > 0) MPI_Send(FIELD+nw, 2, MPI_FLOAT_TYPE, PROCESS_RANK-1, nw, MPI_COMM_WORLD); 
	}
	 
#endif
}


void load_info(FILE* fid)
{
    OMEGA0     = 2*M_PI*LIGHT_VELOCITY/load_namedfloat(fid, "LAMBDA_V");

    TMIN      = load_namedfloat(fid, "T_MIN", true, nan("0")); TMIN = load_namedfloat(fid, "T_MIN_OMEGA", true, TMIN*OMEGA0)/OMEGA0; 
    TMAX      = load_namedfloat(fid, "T_MAX", true, nan("0")); TMAX = load_namedfloat(fid, "T_MAX_OMEGA", true, TMAX*OMEGA0)/OMEGA0;	
    
    N_T       = load_namedint   (fid, "N_T");
    TSTEP     = (TMAX - TMIN)/N_T;
    OMEGA         = (float_type*)malloc_ch(sizeof(float_type)*N_T);
    WAVENUMBER    = (f_complex*) malloc_ch(sizeof(f_complex) *N_T);
	HO_DISPERSION = (float_type*)malloc_ch(sizeof(float_type)*N_T);
     
    char dispfile[300];
	load_namedstringn(fid, "DISPERSION_FILE",dispfile, 300);
	
    calculate_omega(N_T,TMIN,TMAX,OMEGA0,OMEGA);
    calculate_wavenumber(dispfile, N_T,  OMEGA,   WAVENUMBER);
    calculate_wavenumber(dispfile,   1, &OMEGA0, &WAVENUMBER0); 

	{FILE* tf = fopen("wavenum.bin","wb"); fwrite(WAVENUMBER, sizeof(f_complex), N_T, tf); fclose(tf);}
	{FILE* tf = fopen("omega.bin","wb"); fwrite(OMEGA, sizeof(float_type), N_T, tf); fclose(tf);}

	GROUP_VELOCITY = calculate_groupvelocity(N_T, OMEGA, (f_complex*)WAVENUMBER, OMEGA0);

	calculate_HOdispersion();
	if (OMEGA0 == OMEGA[0]) OMEGA_MAX = min(OMEGA_MAX, (1-ABSORBTION_LAYER_WIDTH)*OMEGA[N_T/2-1]);
	else OMEGA_MAX = min(OMEGA_MAX, (1-ABSORBTION_LAYER_WIDTH)*OMEGA[N_T-1]); 

#ifdef NO_FIR
	OMEGA_MIN = max(OMEGA_MIN, 1e15);
#endif

	MAX_TOLERANCE = load_namedfloat(fid, "MAX_TOLERANCE", true, MAX_TOLERANCE_DEFAULT);
	RESIZE_MAX_TOLERANCE = load_namedfloat(fid, "RESIZE_MAX_TOLERANCE", true, DEFAULT_RESIZE_MAX_TOLERANCE);
	RESIZE_MIN_TOLERANCE = load_namedfloat(fid, "RESIZE_MIN_TOLERANCE", true, DEFAULT_RESIZE_MIN_TOLERANCE);
	RESIZE_FREQUENCY = load_namedint(fid, "RESIZE_FREQUENCY", true, DEFAULT_RESIZE_FREQUENCY);


#if   TRANSVERSE_DIMENSIONS == 0
	N_X	    = 1;
	XMIN        = 0;
	XMAX        = 0;

	N_Y	    = 1;
	YMIN	    = 0;
	YMAX	    = 0;
	
#elif TRANSVERSE_DIMENSIONS == 1
	N_X			= load_namedint(fid, "N_R");
	XMIN        = 0;
	XMAX        = load_namedfloat(fid, "R_MAX");
	N_Y			= 1;
	YMIN		= 0;
	YMAX		= 0;
#elif TRANSVERSE_DIMENSIONS == 2 
	N_X         = load_namedint(fid, "N_X");
	XMIN        = load_namedfloat(fid, "X_MIN");
	XMAX        = load_namedfloat(fid, "X_MAX");
	XSTEP       = (XMAX - XMIN)/(N_X); 
	
	N_Y         = load_namedint(fid, "N_Y");
	YMIN        = load_namedfloat(fid, "Y_MIN");
	YMAX        = load_namedfloat(fid, "Y_MAX");
	YSTEP       = (YMAX - YMIN)/(N_Y);
#endif 
	
    N_Z          =  load_namedint(fid, "N_Z");
    float_type zmin  =  load_namedfloat(fid, "Z_MIN", true, 0);
    float_type zmax  =  load_namedfloat(fid, "Z_MAX", true, 0);

    char znettype[50];
    load_namedstringn(fid, "ZNET_TYPE", znettype, 50);
    ZNET  = new float_type[N_Z];
    create_net(zmin, zmax, N_Z, znettype, ZNET);

    CURRENT_Z = ZNET[0];		
	
    APERTURE_N = load_namedfloat(fid, "APERTURE_N", true, 0);
    APERTURE_Z = (float_type*)malloc_ch(sizeof(float_type)*APERTURE_N); 
    APERTURE_R = (float_type*)malloc_ch(sizeof(float_type)*APERTURE_N); 
    for (int i=0; i<APERTURE_N; i++) { APERTURE_Z[i]=load_namednumfloat(fid, "APERTURE_Z",APERTURE_N-i-1); APERTURE_R[i]=load_namednumfloat(fid, "APERTURE_R", APERTURE_N-i-1);}

    AMBIENT_CARRIER_DENSITY = load_namedfloat(fid, "AMBIENT_CARRIERS",true,0)*load_namedfloat(fid, "NEUTRAL_DENSITY",true,0);
}

void  create_mystartcondition(FILE* fid)
{
    #ifndef _SILENCE
	 if (ISMASTER) printf("\n: Creating starting condition...");
    #endif
    char pulseshape[200] = ""; 
    load_namedstringn(fid, "PULSE_SHAPE", pulseshape, 200);
   
    if (ISMASTER)
    {
     printf("\npulse shape = %s", pulseshape);
     fflush(stdout);
    }
    
    if      (strncmp(SHAPE_GG, pulseshape,strlen(SHAPE_GG)) == 0) 									create_gg(fid);
    else if (isdigit(pulseshape[0]) && (strncmp(SHAPE_BEAMSGG, pulseshape+1, strlen(SHAPE_BEAMSGG))==0))        			create_beamsgg(fid);
    else if (isdigit(pulseshape[0]) && (strncmp(SHAPE_BEAMSPROBEGG, pulseshape+1, strlen(SHAPE_BEAMSPROBEGG))==0))			create_beamsprobegg(fid);
    else if (isdigit(pulseshape[0]) && (strncmp(SHAPE_BEAMSPROBEEXTGG, pulseshape+1, strlen(SHAPE_BEAMSPROBEGG))==0))	                create_beamsprobeextgg(fid);
    else if (isdigit(pulseshape[0]) && (strncmp(SHAPE_BEAMSPROBEVGG, pulseshape+1, strlen(SHAPE_BEAMSPROBEVGG))==0))    		create_beamsprobevgg(fid);
    else if (isdigit(pulseshape[0]) && (strncmp(SHAPE_BEAMS_DELAYEDPROBE_GG, pulseshape+1, strlen(SHAPE_BEAMS_DELAYEDPROBE_GG))==0)) 	create_beams_delayedprobegg(fid);
    else if (isdigit(pulseshape[0]) && (strncmp(SHAPE_BEAMS_DELAYEDPROBE_RANDPHASE, pulseshape+1, strlen(SHAPE_BEAMS_DELAYEDPROBE_RANDPHASE))==0)) 	create_beams_delayedprobegg_randphase(fid);
	else if (strncmp(SHAPE_CUSTOMSPECTRUMG, pulseshape, strlen(SHAPE_CUSTOMSPECTRUMG))==0)  		create_customspectrumg(fid);
	else throw "Unknown pulse shape specified!";
}



void load_nonlindata(FILE* fid)
{
	float_type P = load_namedfloat(fid, "PRESSURE", true, 1); 
	NONLIN_REFRINDEX 	= load_namedfloat(fid, "NONLIN_REFRINDEX", true, 0)*P; 
	NONLIN_REFRINDEX4    = load_namedfloat(fid, "NONLIN_REFRINDEX4", true, 0)*P;


	RAMAN_FRACTION   	= load_namedfloat(fid, "RAMAN_FRACTION");
	TAU_RAMAN 	 	= load_namedfloat(fid, "TAU_RAMAN");
	OMEGA_RAMAN      	= load_namedfloat(fid, "OMEGA_RAMAN");
#ifdef THIRD_HARMONICS
	TH_FACTOR		= load_namedfloat(fid, "TH_FACTOR", true, 1);
#endif


	NEUTRAL_DENSITY    	= load_namedfloat(fid, "NEUTRAL_DENSITY")*P;
	RECOMBINATION_TAU  	= load_namedfloat(fid, "RECOMBINATION_TAU"); 
	COLLISION_TAU        	= load_namedfloat(fid, "COLLISION_TAU")/P;
#ifndef MULTI_LEVEL_IONIZATION
	IONIZATION_POTENTIAL 	= load_namedfloat(fid, "IONIZATION_POTENTIAL"); IONIZATION_POTENTIAL *= ELECTRON_CHARGE;
#else
	IONIZATION_LEVEL_N   = load_namedfloat(fid, "IONIZATION_LEVEL_N", true, 1); 
	IONIZATION_POTENTIALS =  (float_type*)malloc_ch(IONIZATION_LEVEL_N*sizeof(float_type)); 
	IONIZATION_RATES_BUFFER = (float_type*)malloc_ch(2*IONIZATION_LEVEL_N*sizeof(float_type)); 
	ION_DENSITIES_BUFFER    = (float_type*)malloc_ch(IONIZATION_LEVEL_N*sizeof(float_type)); 
	IONIZATION_RATE_LN      = (float_type*)malloc_ch(IONIZATION_N*IONIZATION_LEVEL_N*sizeof(float_type));

	N2_IONFACTOR =  (float_type*)malloc_ch(IONIZATION_LEVEL_N*sizeof(float_type)); 

	for (int i=0; i<IONIZATION_LEVEL_N; i++) IONIZATION_POTENTIALS[i] = load_namednumfloat(fid, "IONIZATION_POTENTIAL",i+1)*ELECTRON_CHARGE;
	N2_IONFACTOR[IONIZATION_LEVEL_N-1]=1.0; 
	for (int i=0; i<IONIZATION_LEVEL_N-1; i++) 
  	{
 	 N2_IONFACTOR[i] = load_namednumfloat(fid, "N2_IONFACTOR", i+1, true, 0); 
	 N2_IONFACTOR[IONIZATION_LEVEL_N-1-i] -= N2_IONFACTOR[i]; 
	}

	IONIZATION_POTENTIAL = IONIZATION_POTENTIALS[0]; 	
	
#endif	

#ifdef MULTIPHOTON_IONIZATION
	BETA_MPI_LN		= load_namedfloat(fid, "MPI_CROSSSECTION_LN");
	K_MPI                   = (int)(ceil(IONIZATION_POTENTIAL/(PLANCK_CONSTANT_REDUCED*OMEGA0)));
#endif
#ifdef TUNNEL_IONIZATION
	TUNNELING_FIELD = load_namedfloat(fid, "TUNNELING_FIELD");
#endif

#ifdef BLOCH_RESPONSE
	BLOCH_N = load_namedint(fid, "BLOCH_LEVELS"); 
	BLOCH_TEMP = load_namedfloat(fid, "BLOCH_TEMPERATURE", true, 0); 
	BLOCH_ENERGY = (float_type*)malloc_ch(BLOCH_N*sizeof(float_type)); 
	BLOCH_DIPOLE = (float_type*)malloc_ch(BLOCH_N*BLOCH_N*sizeof(float_type)); 
	BLOCH_G1 = (float_type*)malloc_ch(BLOCH_N*(BLOCH_N)*sizeof(float_type)); 
	BLOCH_G2 = (float_type*)malloc_ch(BLOCH_N*(BLOCH_N)*sizeof(float_type)); 

	float_type phase_v = OMEGA0/real(WAVENUMBER0);
 	float_type C = sqrt(2*phase_v/GROUP_VELOCITY*phase_v*VACUUM_PERMEABILITY);
	float_type DIPOLE_FACTOR = -sqrt(3.0/2.0/ELECTRON_MASS/PLANCK_CONSTANT_REDUCED)*FIELD_DENOM*C*ELECTRON_CHARGE;

	for (int i=0; i<BLOCH_N; i++) BLOCH_ENERGY[i]=load_namednumfloat(fid, "BLOCH_ENERGY", i)*ELECTRON_CHARGE/PLANCK_CONSTANT_REDUCED; //convert energy from eV to s^-1
	for (int i=0; i<BLOCH_N; i++)
	{
	  for (int j=0; j<BLOCH_N; j++) BLOCH_DIPOLE[i+BLOCH_N*j]=0;
	  BLOCH_G1[i + BLOCH_N*i]=0; 
	  BLOCH_G2[i + BLOCH_N*i]=0;
	  for (int j=0; j<i; j++) 
	  {
		
		int nij = i+BLOCH_N*j, nji = j+BLOCH_N*i;
		float_type tau = load_named2numfloat(fid, "BLOCH_T1", j, i, true, 0);
		if (tau > 0) BLOCH_G1[nji]=1.0/tau; else BLOCH_G1[nji]=0; 
		BLOCH_G1[nij]=0;

		tau = load_named2numfloat(fid, "BLOCH_T2", j, i, true, 0);
		if (tau > 0) BLOCH_G2[nij]=1.0/tau; else BLOCH_G2[nij]=0; 
		BLOCH_G2[nji]=BLOCH_G2[nij]; 
		
		BLOCH_DIPOLE[nij] = load_named2numfloat(fid, "BLOCH_DIPOLE", j, i)*DIPOLE_FACTOR/sqrt(BLOCH_ENERGY[i]-BLOCH_ENERGY[j]);
		BLOCH_DIPOLE[nji] = BLOCH_DIPOLE[nij];    
	  }
	}
#endif
}

float_type load_namedfloat(FILE* fid,const char* name, bool defvalue_present, float_type defvalue)
{
    //first, look for the parameter in command line.   
    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0 ) return atof(_ARGV[i+1]);
    
    //then, in the input file.
	fseek(fid, 0, SEEK_SET);
	char namebuf[500];
	char valbuf[300];
	while (!feof(fid))
	{

		fscanf(fid, "%s %s",namebuf, valbuf);
		if (strcmp(namebuf, name) == 0) return atof(valbuf);
	}
	if (!defvalue_present)
	{
		printf("\nload_namedparameter: Parameter %s not found.", name);
		throw "load_namedparameter: parameter not found";
	}
	else  return defvalue;
	
}

float_type load_namednumfloat(FILE* fid,const char* name, int num, bool defvalue_present, float_type defvalue)
{
   //first, look for the parameter in command line.   
    char namenum[500]; sprintf(namenum,"%s%d",name,num);
    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],namenum)==0 ) return atof(_ARGV[i+1]);
    //then, in the input file.
    fseek(fid, 0, SEEK_SET);
    char namebuf[500];
    char valbuf[300];
    while (!feof(fid))
    {
     fscanf(fid, "%s %s",namebuf, valbuf);
     if (strcmp(namebuf, namenum) == 0) return atof(valbuf);
    }

    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0 ) return atof(_ARGV[i+1]);
    //then, in the input file.
    fseek(fid, 0, SEEK_SET);
    while (!feof(fid))
    {
     fscanf(fid, "%s %s",namebuf, valbuf);
     if (strcmp(namebuf, name) == 0) return atof(valbuf);
    }
    if (!defvalue_present)
    {
     printf("\nload_namedparameter: Parameter %s not found.", name);
     throw "load_namedparameter: parameter not found";
    }
    else  return defvalue;
}

float_type load_named2numfloat(FILE* fid,const char* name, int num1, int num2, bool defvalue_present, float_type defvalue)
{
   //first, look for the parameter in command line.   
    char namenum[500]; sprintf(namenum,"%s_%d_%d",name,num1, num2);
    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],namenum)==0 ) return atof(_ARGV[i+1]);
    //then, in the input file.
    fseek(fid, 0, SEEK_SET);
    char namebuf[500];
    char valbuf[300];
    while (!feof(fid))
    {
     fscanf(fid, "%s %s",namebuf, valbuf);
     if (strcmp(namebuf, namenum) == 0) return atof(valbuf);
    }

    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0 ) return atof(_ARGV[i+1]);
    //then, in the input file.
    fseek(fid, 0, SEEK_SET);
    while (!feof(fid))
    {
     fscanf(fid, "%s %s",namebuf, valbuf);
     if (strcmp(namebuf, name) == 0) return atof(valbuf);
    }
    if (!defvalue_present)
    {
     printf("\nload_namedparameter: Parameter %s not found.", name);
     throw "load_namedparameter: parameter not found";
    }
    else  return defvalue;
}




int load_namedint(FILE* fid,const char* name, bool defvalue_present, int defvalue)
{
    //first, look for the parameter in command line.   
    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0) return atoi(_ARGV[i+1]);
    
    //then, in the input file.

	fseek(fid, 0, SEEK_SET);
	char namebuf[500];
	char valbuf[300];
	while (!feof(fid))
	{

		fscanf(fid, "%s %s",namebuf, valbuf);
		if (strncmp(namebuf, name, 500) == 0) return atoi(valbuf);
	}
	if (!defvalue_present)
	{
		printf("\nload_namedint: Parameter %s not found.", name);
		throw "load_namedint: parameter not found";
	}
	else return defvalue;
}

int load_namednumint(FILE* fid,const char* name_, int num, bool defvalue_present, int defvalue)
{
    //first, look for the parameter in command line.   
	char name[500]; sprintf(name,"%s%d",name_,num);
    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0) return atoi(_ARGV[i+1]);
    
    //then, in the input file.

	fseek(fid, 0, SEEK_SET);
	char namebuf[500];
	char valbuf[300];
	while (!feof(fid))
	{

		fscanf(fid, "%s %s",namebuf, valbuf);
		if (strncmp(namebuf, name, 500) == 0) return atoi(valbuf);
	}
	if (!defvalue_present)
	{
		printf("\nload_namedint: Parameter %s not found.", name);
		throw "load_namedint: parameter not found";
	}
	else return defvalue;
}

void load_namedstringn(FILE* fid,const char* name, char* output, int N, bool defvalue_present, const char* defvalue)
{
    //first, look for the parameter in command line.   
    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0 ) {strncpy(output, _ARGV[i+1], N); return;}
    
    //then, in the input file.
 
	fseek(fid, 0, SEEK_SET);
	char namebuf[500];
	char valbuf[300];
	while (!feof(fid))
	{
		fscanf(fid, "%s %s",namebuf, valbuf);
		if (strncmp(namebuf, name, 500) == 0) {strncpy(output, valbuf, N); return;}
	}
	if (!defvalue_present)
	{
		printf("\nload_nameddouble: Parameter %s not found.", name);
		throw "load_namedfloat: parameter not found";
	}
	else strncpy(output, defvalue, N); 
}


void load_namednumstringn(FILE* fid,const char* name_, char* output, int num, int N, bool defvalue_present, const char* defvalue)
{
    //first, look for the parameter in command line.   
	char name[500]; sprintf(name,"%s%d",name_,num);
    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0 ) {strncpy(output, _ARGV[i+1], N); return;}
    
    //then, in the input file.
 
	fseek(fid, 0, SEEK_SET);
	char namebuf[500];
	char valbuf[300];
	while (!feof(fid))
	{

		fscanf(fid, "%s %s",namebuf, valbuf);
		if (strncmp(namebuf, name, 500) == 0) {strncpy(output, valbuf, N); return;}
	}
	if (!defvalue_present)
	{
		printf("\nload_nameddouble: Parameter %s not found.", name);
		throw "load_nameddouble: parameter not found";
	}
	else strncpy(output, defvalue, N); 
}



void print_variables()
{
    printf(  "\nTIME_START  =%e", TIME_START);

#ifdef _FULL_DUMPS
    printf(  "\nDUMP_ID     =%d", DUMP_ID);
    printf(  "\nDUMP_PREFIX =%s", DUMP_PREFIX);
#endif
    
    printf("\n\nTMIN=%e, TMAX=%e, N_T=%d", TMIN, TMAX, N_T);
#ifdef _UNWRAP_FREQUENCIES
	printf("\n Unwrapping frequencies : on");
#endif

	printf("\n Maximum ZSTEP over nonlinear length: %e", MAX_TOLERANCE);

    printf(  "\nN_X=%d, XMIN=%e, XMAX=%e",N_X, XMIN, XMAX);
    printf(  "\nN_Y=%d, YMIN=%e, YMAX=%e",N_Y, YMIN, YMAX); 
    
    printf("\nN_Z=%d",N_Z);
    printf("\nOMEGA0 = %e, WAVENUMBER0 = %e + %ei", (double)OMEGA0, (double)real(WAVENUMBER0), (double)imag(WAVENUMBER0));
    printf("\nGROUP_VELOCITY=%e", (double)GROUP_VELOCITY);

	printf("\n\n INTENSITY_DENOM = %e", INTENSITY_DENOM);
	printf("\n\n IONRATE_DENOM = %e", IONRATE_DENOM);

    printf("\n\nNONLIN_REFRINDEX  = %e", NONLIN_REFRINDEX/INTENSITY_DENOM);
	printf("\n\nNONLIN_REFRINDEX4 (cm^4/TW^2) = %e", NONLIN_REFRINDEX4/INTENSITY_DENOM/INTENSITY_DENOM*1e16*1e16);
    printf(  "\nRAMAN_FRACTION   = %e", RAMAN_FRACTION);
    printf(  "\nTAU_RAMAN        = %e", TAU_RAMAN);
    printf(  "\nOMEGA_RAMAN      = %e", OMEGA_RAMAN);

#ifdef THIRD_HARMONICS
	printf("\n THIRD_HARMONICS: on; TH_FACTOR = %g", TH_FACTOR);
#else 
	printf("\n THIRD_HARMONICS: off");
#endif

    printf("\n\nNEUTRAL_DENSITY        = %e", NEUTRAL_DENSITY);
    printf(  "\nAVALANCHE_CROSSSECTION = %e", AVALANCHE_CROSSSECTION);
    printf(  "\nRECOMBINATION_TAU      = %e", RECOMBINATION_TAU/IONRATE_DENOM);
    printf(  "\nCOLLISION_TAU          = %e", COLLISION_TAU);
    printf(  "\nIONIZATION_POTENTIAL   = %g", IONIZATION_POTENTIAL*INTENSITY_DENOM/IONRATE_DENOM/ELECTRON_CHARGE);
	printf(  "\nPONDEROMOTIVE_K        = %e", PONDEROMOTIVE_COEFFICIENT/IONRATE_DENOM);
#ifdef MULTIPHOTON_IONIZATION
    printf(  "\nBETA_MPI_LN            = %e", BETA_MPI_LN);
    printf(  "\nK_MPI                  = %d", K_MPI);
#endif
#ifdef TUNNEL_IONIZATION
	printf(  "\nTUNNELING_FIELD   = %e", TUNNELING_FIELD);
#endif

#ifdef MULTI_LEVEL_IONIZATION
	for (int i=0; i<IONIZATION_LEVEL_N; i++) 
	printf(   "\nIONIZATION_POTENTIAL%d  = %g, n2 factor = %g", i+1, IONIZATION_POTENTIALS[i]*INTENSITY_DENOM/IONRATE_DENOM/ELECTRON_CHARGE, N2_IONFACTOR[i]);
#endif
#ifdef YUDIN_IVANOV_CORRECTION
	printf("\n Yudin-Ivanov correction for instantaneous ionization rate is on.");
#endif
	fflush(stdout);
}

void create_net(float_type Xmin, float_type Xmax, int N, char* nettype, float_type* net)
{
	if (nettype == NULL || nettype[0]== 'e')
	{
		//Equi-step net
		float_type step = (Xmax - Xmin)/(N-1);
		for (int i=0; i<N; i++)	net[i] = Xmin + i*step;
	}
	else if (nettype[0]=='p')
	{
		//Odd power net
		int pw = atoi(nettype+1); 
        float_type n0 = N/(1-oddroot(Xmax/Xmin,pw));
        float_type A = -Xmin/oddpow(n0,pw);
        for (int i=0; i<N; i++) net[i] = A*oddpow(i-n0,pw);
	}

	else throw "Unrecognized net type!";
}


void calculate_omega(int Nt, float_type tmin, float_type tmax, float_type omega0, float_type* omega)
{
	float_type wstep = 2*M_PI/(tmax-tmin);
	float_type omega_hw = wstep*Nt/2;

#ifdef _UNWRAP_FREQUENCIES
	if (omega0 <= omega_hw) for (int i=0;     i<Nt; i++) omega[i] = wstep*(i+1);
	else
#endif 
	{
		for (int i=0;     i<Nt/2; i++) omega[i] = omega0 + wstep*i;
		for (int i=Nt/2;  i<Nt;   i++) omega[i] = omega0 - wstep*(Nt-i);
	}
}


void calculate_wavenumber(char* filename, int N, float_type* omega, f_complex* wavenum)
{
	int filetype = 0;
	FILE* fid = NULL;
#ifdef _DEBUG
	printf("\ncalculate_wavenumber : loading dispersion information from file %s...", filename);
#endif
	if (filename[0] == 0 || !strcmp(filename,"NO"))
	{
		filetype = FILETYPE_NOFILE;
	}
	else
	{
		fid = fopen(filename, "rb"); 
		if (!fid) throw "Unable to open file with refractive index information!";
		fread(&filetype, sizeof(int), 1, fid);
	}
	f_complex* refr_index = (f_complex*)malloc_ch(sizeof(f_complex)*N);
	switch (filetype)
	{
		case FILETYPE_NOFILE           : for (int i=0; i<N; i++) refr_index[i] = 1.0;                         break;
		case FILETYPE_SELLMEIER_LAMBDA : load_refrindex_sellmeier_lambda(fid, N, omega, refr_index); break;
		case FILETYPE_SELLMEIER_OMEGA  : load_refrindex_sellmeier_omega (fid, N, omega, refr_index); break;
		case FILETYPE_RAW              : load_refrindex_raw             (fid, N, omega, refr_index); break;
		case FILETYPE_DISPROW          : {load_dispersion_row(fid, N, omega, wavenum); free(refr_index); fclose(fid); return;} 
		default: throw "Unknown type of file with refractive index information!";
	}
	float_type omega_max=omega[0]; for (int i=0; i<N; i++) omega_max = max(omega_max,omega[i]);

	for (int i=0; i<N; i++) 
	{	
	  wavenum[i]  = conj(refr_index[i])*omega[i]/LIGHT_VELOCITY;
         f_complex j = f_complex(0,1);
         float_type domega = (omega[i]/omega_max - 1 - ABSORBTION_LAYER_WIDTH)/ABSORBTION_LAYER_WIDTH;
         if (domega > 0) wavenum[i] += -j*real(wavenum[i])*(exp_p(ABSORBTION_LAYER_BETA)*domega*domega);
    	}
	free(refr_index);
	if (filetype != FILETYPE_NOFILE) fclose(fid);
}



void load_refrindex_sellmeier_omega (FILE* fid, int N, float_type* cfreq, f_complex* refrindex)
{
	int koefN = 0;
	float_type* sB = NULL;
	float_type* somega = NULL;
	fread(&koefN, sizeof(int), 1, fid);

	sB      = new(float_type[koefN]);
	somega  = new(float_type[koefN]);
	fread(somega,  sizeof(float_type), koefN, fid);
	fread(sB,      sizeof(float_type), koefN, fid);
		
	for (int i=0; i<N; i++)
	{
		f_complex n2 = 1.0;
		for (int j=0; j<koefN; j++) n2 += sB[j]*(somega[j]*somega[j])/(somega[j]*somega[j] - cfreq[i]*cfreq[i]);
		refrindex[i] = sqrt(n2);
	}
	delete sB;
	delete somega;
}

void load_refrindex_sellmeier_lambda(FILE* fid, int N, float_type* cfreq, f_complex* refrindex)
{
	int koefN = 0;
	float_type* sB	= NULL;
	float_type* slambda = NULL;
	float_type* somega  = NULL;
	fread(&koefN, sizeof(int), 1, fid);

	sB      = new(float_type[koefN]);
	somega  = new(float_type[koefN]);
	slambda = new(float_type[koefN]);

	fread(slambda,  sizeof(float_type), koefN, fid);
	fread(sB,      sizeof(float_type), koefN, fid);
		
	for (int j=0; j<koefN; j++) somega[j] = 2*M_PI*LIGHT_VELOCITY/slambda[j];

	for (int i=0; i<N; i++)
	{
		f_complex n2 = 1.0;
		for (int j=0; j<koefN; j++) n2 += sB[j]*(somega[j]*somega[j])/(somega[j]*somega[j] - cfreq[i]*cfreq[i]);
		refrindex[i] = sqrt(n2);	
	}
	delete sB;
	delete slambda;
	delete somega;
}

void load_refrindex_raw(FILE* fid, int N, float_type* cfreq, f_complex* refrindex)
{
#ifdef _DEBUG
	printf("\n void load_refrindex_raw(FILE*, int, float_type*, float_type*)");
	printf("\n N = %d", N);
#endif
	float_type P = load_namedfloat(SC_FID, "PRESSURE",true,1);
	int pointsN = 0;
	double*      pomega = NULL;
	double*      pn     = NULL;
	fread(&pointsN,sizeof(int),1, fid);

	pomega = (double*)malloc_ch(sizeof(double)*pointsN);
	pn     = (double*)malloc_ch(2*sizeof(double)*pointsN);
	
	fread(pomega,sizeof(double),   pointsN, fid);
	fread(pn    ,sizeof(double), 2*pointsN, fid);


	OMEGA_MAX = pomega[pointsN-1]; OMEGA_MIN = pomega[0];
	for (int j=0; j<pointsN; j++) {OMEGA_MAX = max(OMEGA_MAX,pomega[j]); OMEGA_MIN = min(OMEGA_MIN, pomega[j]);}
	for (int i=0; i<N; i++)
	{
		for (int j=1; j<pointsN; j++)
		{
			if ((pomega[j-1]-cfreq[i])*(pomega[j]-cfreq[i])<=0) 
				{
					double ren_ = P*(pn[2*j-2]-1)+1, imn_ = P*pn[2*j-1];
					double ren  = P*(pn[2*j]  -1)+1, imn  = P*pn[2*j+1];
					double renc = ren_ + (ren-ren_)*(fabs(cfreq[i])-pomega[j-1])/(pomega[j]-pomega[j-1]); 
					double imnc = imn_ + (imn-imn_)*(fabs(cfreq[i])-pomega[j-1])/(pomega[j]-pomega[j-1]);

					refrindex[i] = f_complex(renc, imnc);
					break;
				}
			if (j==(pointsN-1))
			{
				refrindex[i]=1;
			}
		}
	}

    free(pomega);
	free(pn);
}


float_type calculate_groupvelocity(int Nt, float_type* omega, f_complex* wavenum, float_type omega0)
{
	float_type wstep = fabs(omega[1]-omega[0]);

	for (int i=0; i<Nt-1; i++) if (fabs(omega0-omega[i])<wstep) return (omega[i]-omega[i+1])/real(wavenum[i]-wavenum[i+1]);
	throw "calculate_groupvelocity: invalid input: omega0 is not inside the omega net";	

}


void init_zstep_kerr()
{
 float_type maxI = 0;
 for (int ny=0; ny<MY_NY; ny++)
 for (int nx=0; nx<MY_NX;   nx++)
 {
  for (int nt=0; nt<N_T; nt++)
  {
   int ofs =nt+N_T*(nx+MY_NX*ny);
   maxI = max(maxI, abs2(BIGBUFFER1[ofs]));
  }
 }
 float_type maxI_ = maxI;
 MPI_Allreduce(&maxI_,  &maxI,  1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  

 ZSTEP = ZNET[1]-ZNET[0];
 float_type kNL = OMEGA0/LIGHT_VELOCITY*NONLIN_REFRINDEX*maxI;
 if (ZSTEP*kNL > MAX_TOLERANCE) ZSTEP = MAX_TOLERANCE/kNL;
 if (ISMASTER) {printf("\n Initializing ZSTEP according to Kerr nonlinearity. maxI=%e, kNL=%e, ZSTEP=%e", maxI, kNL, ZSTEP); fflush(stdout);}
}


#ifdef _EXTMATH_SINGLE

void load_refrindex_raw(FILE* fid, int N, float_type* cfreq, d_complex* refrindex)
{
#ifdef _DEBUG
	printf("\n void load_refrindex_raw(FILE*, int, float_type*, float_type*)");
	printf("\n N = %d", N);
#endif
	float_type P = load_namedfloat(SC_FID, "PRESSURE",true,1);
	int pointsN = 0;
	double*      pomega = NULL;
	double*      pn     = NULL;
	fread(&pointsN,sizeof(int),1, fid);

	pomega = (double*)malloc_ch(sizeof(double)*pointsN);
	pn     = (double*)malloc_ch(2*sizeof(double)*pointsN);
	
	fread(pomega,sizeof(double),   pointsN, fid);
	fread(pn    ,sizeof(double), 2*pointsN, fid);


	OMEGA_MAX = pomega[pointsN-1]; OMEGA_MIN = pomega[0];
	for (int j=0; j<pointsN; j++) {OMEGA_MAX = max(OMEGA_MAX,pomega[j]); OMEGA_MIN = min(OMEGA_MIN, pomega[j]);}
	for (int i=0; i<N; i++)
	{
		for (int j=1; j<pointsN; j++)
		{
			if ((pomega[j-1]-cfreq[i])*(pomega[j]-cfreq[i])<=0) 
				{
					double ren_ = P*(pn[2*j-2]-1)+1, imn_ = P*pn[2*j-1];
					double ren  = P*(pn[2*j]  -1)+1, imn  = P*pn[2*j+1];
					double renc = ren_ + (ren-ren_)*(fabs(cfreq[i])-pomega[j-1])/(pomega[j]-pomega[j-1]); 
					double imnc = imn_ + (imn-imn_)*(fabs(cfreq[i])-pomega[j-1])/(pomega[j]-pomega[j-1]);

					refrindex[i] = d_complex(renc, imnc);
					break;
				}
			if (j==(pointsN-1))
			{
				refrindex[i]=1;
			}
		}
	}

    free(pomega);
	free(pn);
}

#endif

void calculate_HOdispersion()
{
	char dispfilename		[300]; 
	load_namedstringn(SC_FID, "DISPERSION_FILE", dispfilename, 300, true, "NO");
	if (strncmp(dispfilename, "NO", 300) == 0) {for (size_t nt=0; nt<N_T; nt++) HO_DISPERSION[nt] = 0.0; return;}

	FILE* dispfid = fopen(dispfilename, "rb"); if (!dispfid) throw "Error opening refractive index file!";
	int type=0;
	fread(&type, sizeof(int), 1, dispfid);
	if (type == FILETYPE_RAW) 
	{
	 d_complex* refr_index = (d_complex*)malloc(N_T*sizeof(d_complex)); 
	 load_refrindex_raw(dispfid, N_T, OMEGA, refr_index);
	 double wstep = fabs(OMEGA[1]-OMEGA[0]);
	 double n0=0, dndw = 0; 
	 if (OMEGA0 != OMEGA[0])
	 {
	   for (int i=0; i<N_T-1; i++) 
	      if (fabs(OMEGA0-OMEGA[i])<wstep) 
		  {
			 dndw = real(refr_index[i+1]-refr_index[i])/((double)OMEGA[i+1]-(double)OMEGA[i]);
			 n0 = real(refr_index[i]) + dndw*(OMEGA0-OMEGA[i]);
			 break;
		  }
	 }
	 else
	 {
		  n0 = real(refr_index[0]);
		  dndw = real(refr_index[1]-refr_index[N_T-1])/((double)OMEGA[1]-(double)OMEGA[N_T-1]);
	 }

	 double ngr = n0 + OMEGA0*dndw;
			 	 
	 for (size_t nw=0; nw < N_T; nw++) 
	 { 
		HO_DISPERSION[nw] = (float_type)((OMEGA[nw])/dLIGHT_VELOCITY*(real(refr_index[nw]) - ngr));
	 }

	 GROUP_VELOCITY = LIGHT_VELOCITY/ngr;

	}
	else if (type == FILETYPE_DISPROW)
	{
		float_type P = load_namedfloat(SC_FID, "PRESSURE", true, 1);
		int N; fread(&N, sizeof(int), 1, dispfid); 
		double omega0; fread(&omega0, sizeof(double), 1, dispfid);
		if (abs(omega0-OMEGA0)/OMEGA0 > 1e-3) throw "Incorrect dispersion expansion! Center expansion frequency must coincide with pulse frequency";
		double* dispk = (double*)malloc_ch(sizeof(double)*N);  fread(dispk, sizeof(double), N, dispfid); fclose(dispfid);
		WAVENUMBER0 = (f_complex)((dispk[0]-((double)OMEGA0)/(double)LIGHT_VELOCITY)*P+(double)OMEGA0/(double)LIGHT_VELOCITY);
		GROUP_VELOCITY = 1.0/((dispk[1]-1.0/(double)LIGHT_VELOCITY)*P + 1.0/(double)LIGHT_VELOCITY);
		for (size_t nw=0; nw<N_T; nw++)
		{
			HO_DISPERSION[nw] = 0; 
			for (int i=2; i<N; i++) HO_DISPERSION[nw]+=P*dispk[i]*pow((double)OMEGA[nw]-(double)OMEGA0,i);
			WAVENUMBER[nw] = WAVENUMBER0 + (f_complex)((double)OMEGA[nw]-(double)OMEGA0)/GROUP_VELOCITY + HO_DISPERSION[nw];
		}
		free(dispk);
	}
	else throw "HO-dispersion calculation from refractive index data of other type than raw is not supported!";
	for (size_t nw=0; nw<N_T; nw++)
	{
		float_type w = OMEGA[nw]; 
		float_type f = lfgaussfilter(w)*hfgaussfilter(w); 
		HO_DISPERSION[nw]*=f;
	}
}


void load_spectrum_fromfile(f_complex* output, float_type* omega, int Nt, const char* spectrumfilename, const char* spectrumfiletype)
{

	
	f_complex j = f_complex(0,1);

	int N=0; 
	float_type *somega, *sI, *sphi; 

	if (strncmp(spectrumfiletype, "text_lambdanm",128) == 0)
	{
		FILE* sfid = fopen(spectrumfilename, "rt"); if (!sfid) throw "Unable to open input spectrum file!";
		char buf[512];
		fgets(buf, 512, sfid); 
		float_type somega_, sI_, sphi_; 

#ifdef _EXTMATH_SINGLE
		char format3f[] = "%f %f %f"; char format2f[] = "%f %f"; 
#else
		char format3f[] = "%lf %lf %lf"; char format2f[] = "%lf %lf"; 
#endif

		int ncols = sscanf(buf, format3f, &somega_, &sI_, &sphi_);  
		if (ncols < 2 || ncols > 3) throw "Incorrect spectrum file for text_lambdanm type!";
		while (!feof(sfid)) {fgets(buf, 512, sfid); N++;} rewind(sfid);

		somega = (float_type*)malloc(sizeof(float_type)*N);
		sI     = (float_type*)malloc(sizeof(float_type)*N);
		sphi   = (float_type*)malloc(sizeof(float_type)*N);

		for (int i=0; i<N; i++)
		{
			if (ncols == 2) {fscanf(sfid, format2f, somega+i, sI+i); sphi[i]=0;}
			else		     fscanf(sfid, format3f, somega+i, sI+i, sphi+i); 
			somega[i] = 2.0*M_PI*LIGHT_VELOCITY/somega[i]/1e-9;
			if (sI[i] < 0) sI[i]=0; 
		}
		fclose(sfid);
	}
	else throw "Unsupported spectrum type!";

	// linear interpolation part
	float_type sE = 0;  //this variable will contain spectrum energy
	for (int nw=0; nw<Nt; nw++)
	{
		int i=0; 
		float_type I, phi; 
		for (; i<N-1; i++) 	
		{
			if ((somega[i] == omega[nw])) 
			{
				I = sI[i];
				phi = sphi[i]; 
				output[nw] = ((float_type)sqrt(I))*exp(-j*phi);
				break;
			}
			if ((somega[i]-omega[nw])*(somega[i+1]-omega[nw]) < 0)
			{
				I   = sI[i]   + (  sI[i+1]  -sI[i])*(omega[nw]-somega[i])/(somega[i+1]-somega[i]);
				phi = sphi[i] + (sphi[i+1]-sphi[i])*(omega[nw]-somega[i])/(somega[i+1]-somega[i]);
				output[nw] = ((float_type)sqrt(I))*exp(-j*phi);
				sE += I;
				break;
			}
		}
		if (i==N-1) output[nw]=0.0;;
	}

	sE = sqrt(sE/Nt);
	for (int nw=0; nw<Nt; nw++) output[nw] /= sE;
	free(somega); free(sI); free(sphi);
}

void load_dispersion_row(FILE* fid, int N, float_type* omega, f_complex* wavenum)
{
	fseek(fid, 0, SEEK_SET);
	int type = 0; fread(&type, sizeof(int), 1, fid);
	if (type != FILETYPE_DISPROW) throw "load_dispersion_row() : Incorrect file type!";
	int Nk = 0; fread(&Nk, sizeof(int), 1, fid);
	double omega0 = 0; fread(&omega0, sizeof(double), 1, fid);
	double* K; K=(double*)malloc_ch(sizeof(double)*Nk);
	fread(K, sizeof(double), Nk, fid);
	for (int nw=0; nw<N; nw++) 
	{
		wavenum[nw] = K[0]; 
		for (int i=1; i<Nk; i++) wavenum[nw]+=K[i]*pow((omega[nw]-omega0), i);
	}

	free(K);
}


void load_kerr_th_profiles(FILE* scfid)
{
  char fname[300]; 
  load_namedstringn(scfid, "KERR_PROFILE", fname, 300, true, "NO"); 
  float_type P = load_namedfloat(scfid, "PRESSURE", true, 1); 
  
#ifndef NO_SHOCK
  for (int nw=0; nw<N_T; nw++) KERR_PROFILE[nw]=NONLIN_REFRINDEX*OMEGA[nw]/LIGHT_VELOCITY/real(WAVENUMBER[nw])*OMEGA[nw]/LIGHT_VELOCITY; 
#else
  for (int nw=0; nw<N_T; nw++) KERR_PROFILE[0]=NONLIN_REFRINDEX*OMEGA[0]/LIGHT_VELOCITY/real(WAVENUMBER[0])*OMEGA[0]/LIGHT_VELOCITY; 
#endif
  
  if (strncmp(fname, "NO", 300))
  {
   if (ISMASTER) printf("\n Loading Kerr profile info from %s .", fname); 
   FILE* fid = fopen(fname, "rb"); 
   int type = 0; fread(&type, 1, sizeof(int), fid); if (type != FILETYPE_KERR_PROFILE) throw "Improper file type!"; 
   int Nt_ = 0; fread(&Nt_, 1, sizeof(int), fid); 
   double *w_=NULL, *n2_=NULL;  w_=(double*)malloc_ch(Nt_*sizeof(double)); n2_=(double*)malloc_ch(Nt_*sizeof(double)); 
   fread(w_, Nt_, sizeof(double), fid); 
   fread(n2_, Nt_, sizeof(double), fid);
   fclose(fid); 

   if (OMEGA_MAX > w_[Nt_-1]) OMEGA_MAX = w_[Nt_-1]; 
   if (OMEGA_MIN < w_[0])     OMEGA_MIN = w_[0];

 /*  gsl_interp* spl          = gsl_interp_alloc(gsl_interp_linear, Nt_);
   gsl_interp_accel* spl_ac =  gsl_interp_accel_alloc(); 
   gsl_interp_init(spl, w_, n2_, Nt_); 
   
   for (int nw=0; nw<Nt_; nw++) {if (OMEGA[nw]>OMEGA_MAX || OMEGA[nw]<OMEGA_MIN) KERR_PROFILE[nw]=0.0; else KERR_PROFILE[nw] *= (float_type)gsl_interp_eval(spl, w_, n2_, OMEGA[nw], spl_ac)*INTENSITY_DENOM/NONLIN_REFRINDEX; } 

   gsl_interp_free(spl); 
   gsl_interp_accel_free(spl_ac); */

   for (int nw=0; nw<N_T; nw++) 
   
    if (OMEGA[nw]>OMEGA_MAX || OMEGA[nw]<OMEGA_MIN) KERR_PROFILE[nw]=0.0; 
    else for (int nw_=1; nw_<=Nt_; nw_++) 
     {
      if (nw_==Nt_) 
          throw "load_kerr_th_profiles(): Interpolation error!";
      if (w_[nw_-1] <= OMEGA[nw] && OMEGA[nw] < w_[nw_])
      {
       double n2m = n2_[nw_-1], n2p = n2_[nw_]; 
       double n2 = n2m + (OMEGA[nw]-w_[nw_-1])/(w_[nw]-w_[nw_-1])*(n2p-n2m); 
       KERR_PROFILE[nw] *= n2*INTENSITY_DENOM*P/NONLIN_REFRINDEX; 		
       break;
      }
     } 
   free(w_); free(n2_); 
  }

  load_namedstringn(scfid, "TH_PROFILE", fname, 300, true, "NO"); 
  for (int nw=0; nw<N_T; nw++) KERR_TH_PROFILE[nw]=KERR_PROFILE[nw]*(1-RAMAN_FRACTION)/3.0*TH_FACTOR;
  
  if (strncmp(fname, "NO", 300))
  {
   if (ISMASTER) printf("\n Loading third-harmonics profile info from %s .", fname); 
   FILE* fid = fopen(fname, "rb"); 
   int type = 0; fread(&type, 1, sizeof(int), fid); if (type != FILETYPE_KERR_TH_PROFILE) throw "Improper file type!"; 
   int Nt_ = 0; fread(&Nt_, 1, sizeof(int), fid); 
   double *w_=NULL, *f_=NULL;  w_=(double*)malloc_ch(Nt_*sizeof(double)); f_=(double*)malloc_ch(Nt_*sizeof(double)); 
   fread( w_, Nt_, sizeof(double), fid); 
   fread(f_, Nt_, sizeof(double), fid);
   fclose(fid); 

   if (OMEGA_MAX > w_[Nt_-1]) OMEGA_MAX = w_[Nt_-1]; 
   if (OMEGA_MIN < w_[0])     OMEGA_MIN = w_[0];

   /*gsl_interp* spl = gsl_interp_alloc(gsl_interp_linear, Nt_);
   gsl_interp_init(spl, w_, n2_, Nt_); 
   gsl_interp_accel* spl_ac =  gsl_interp_accel_alloc(); 
   
   for (int nw=0; nw<Nt_; nw++) {if (OMEGA[nw]>OMEGA_MAX || OMEGA[nw]<OMEGA_MIN) KERR_TH_PROFILE[nw]=0.0; else KERR_TH_PROFILE[nw] *= gsl_interp_eval(spl, w_, n2_, OMEGA[nw], spl_ac); } 
	
   gsl_interp_accel_free(spl_ac); 
   gsl_interp_free(spl);  */
   for (int nw=0; nw<N_T; nw++) 
   
    if (OMEGA[nw]>OMEGA_MAX || OMEGA[nw]<OMEGA_MIN) KERR_TH_PROFILE[nw]=0.0; 
    else for (int nw_=1; nw_<=Nt_; nw_++) 
     {
      if (nw_==Nt_) throw "load_kerr_th_profiles(): Interpolation error!";
      if (w_[nw_-1] <= OMEGA[nw] && OMEGA[nw] < w_[nw_])
      {
       double fm = f_[nw_-1], fp = f_[nw_]; 
       double f = fm + (OMEGA[nw]-w_[nw_-1])/(w_[nw]-w_[nw_-1])*(fp-fm); 
       KERR_TH_PROFILE[nw] *= f;
       break;
      }
     } 
   free(w_); free(f_); 
  }
}


  
