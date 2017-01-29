#include "solver.h"
#include "pulseshapes.h"

void calculate_omega(int Nt, float_type tmin, float_type tmax, float_type omega0, float_type* omega);
void calculate_wavenumber(char* filename, int N, float_type* omega, f_complex* wavenum);

float_type calculate_groupvelocity(int Nt, float_type* omega, f_complex* wavenum, float_type omega0);

void load_refrindex_sellmeier_omega (FILE* fid, int N, float_type* cfreq, f_complex* refrindex);
void load_refrindex_sellmeier_lambda(FILE* fid, int N, float_type* cfreq, f_complex* refrindex);
void load_refrindex_raw             (FILE* fid, int N, float_type* cfreq, f_complex* refrindex);
void load_refrindex_raw_vec			(FILE* fid, int N, float_type* cfreq, f_complex* refrindex);

void create_net(float_type Xmin, float_type Xmax, int N, char* nettype, float_type* net);
float_type integrate_gaussian(float_type xmin, float_type xmax, float_type width);

void init_zstep_kerr();

enum 		{FILETYPE_NOFILE               = 0, \
        	 FILETYPE_SELLMEIER_LAMBDA     = 1, \
			 FILETYPE_SELLMEIER_OMEGA      = 2, \
			 FILETYPE_RAW                  = 3,	\
			 FILETYPE_SELLMEIER_LAMBDA_VEC = 4, \
 			 FILETYPE_SELLMEIER_OMEGA_VEC  = 5, \
			 FILETYPE_RAW_VEC	   = 6};

void process_input(int argc, char** argv)
{
        if (argc < 2 && PROCESS_RANK == 0) throw "\n!!! Error: Too few arguements  !!!";	

	#ifndef _SILENCE
	if (ISMASTER) printf("\n Loading starting info from %s", argv[1]);
	#endif
	FILE* scfid = fopen(argv[1], "r");
	if (scfid == NULL) throw "Error opening file with starting condition.";
	load_info(scfid);
	load_nonlindata(scfid);
	#ifndef _SILENCE
	if (ISMASTER) printf("\n Initializing variables...");
	#endif
	initialize_variables();

    if (ISMASTER) print_variables();	

	DUMPMAN = new dump_manager_class(scfid);	
	
	create_mystartcondition(scfid);	
    
	fclose(scfid);
	// create_mystartcondition() outputs Fourier-transform of the field to FIELD array. 
	// To estimate initial zstep, one needs to perform an inverse transform to BIGBUFFER1
	// Two transforms - one is for first vector component (ordinary wave), another - for second (extraordinary wave)
	fftw(FFT_BWPLAN_T, N_X*MY_NY, (fftwt_complex*)FIELD,     2,2*N_T, (fftwt_complex*)BIGBUFFER1,     2,2*N_T); 
	fftw(FFT_BWPLAN_T, N_X*MY_NY, (fftwt_complex*)(FIELD+1), 2,2*N_T, (fftwt_complex*)(BIGBUFFER1+1), 2,2*N_T);
	fftwt_Nnormalize(2*N_X*MY_NY, BIGBUFFER1);

	ZSTEP = ZNET[1]-ZNET[0];
	init_zstep_kerr();

	memcpy(BIGBUFFER2, FIELD, 2*sizeof(f_complex)*N_T*N_X*MY_NY); 
#ifdef _UNIAXIAL
	//FHT_PLAN->run_many(BIGBUFFER2, 2*N_T, 2*N_T,1);
        fhatha_runmany_cuda(FHT_PLAN, BIGBUFFER2, 2*N_T, 2*N_T, 1); 

#else
	fftwnd_mpi(FFT_FWPLAN_XY, 2*N_T, (fftwt_complex*)BIGBUFFER2, (fftwt_complex*)BIGBUFFER1, FFTW_TRANSPOSED_ORDER);
#endif 
/*	if (ISMASTER)
        {
	 printf("\n Ionization rate:\n ln(I)    ln(W)");
         for (int i=0; i<IONIZATION_N; i++) printf("\n %f       %f", IONIZATION_I_LN[i], IONIZATION_RATE_LN[i]);
         fflush(stdout);
        } */
}


void load_info(FILE* fid)
{
    TMIN      = load_namedfloat(fid, "T_MIN");
    TMAX      = load_namedfloat(fid, "T_MAX");
    
    N_T       = load_namedint   (fid, "N_T");
    TSTEP     = (TMAX - TMIN)/N_T;

	THETA_OA  = load_namedfloat(fid, "THETA_OA");
	PHI_OA    = load_namedfloat(fid, "PHI_OA");

    OMEGA0      = 2*M_PI*LIGHT_VELOCITY/load_namedfloat(fid, "LAMBDA_V");
    OMEGA       = (float_type*)malloc_ch(  sizeof(float_type)*N_T);
    WAVENUMBER  = (f_complex*) malloc_ch(2*sizeof(f_complex)*N_T);
	WAVENUMBER0 = (f_complex*) malloc_ch(2*sizeof(f_complex));
     
    char dispfile[300];
	load_namedstringn(fid, "DISPERSION_FILE",dispfile, 300);
	
    calculate_omega(N_T,TMIN,TMAX,OMEGA0,OMEGA);
    calculate_wavenumber(dispfile, N_T,  OMEGA,   WAVENUMBER);
    calculate_wavenumber(dispfile,   1, &OMEGA0,  WAVENUMBER0); 
	
	OMEGA_MAX = min(OMEGA_MAX, (1-ABSORBTION_LAYER_WIDTH)*OMEGA[N_T/2-1]);

	GROUP_VELOCITY = calculate_groupvelocity(N_T, OMEGA, WAVENUMBER, OMEGA0);

#ifdef _UNIAXIAL
	N_X			= load_namedint(fid, "N_R");
	XMIN        = 0;
	XMAX        = load_namedfloat(fid, "R_MAX");
	N_Y			= 1;
	YMIN		= 0;
	YMAX		= 0;
#else
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
    float_type zmin  =  load_namedfloat(fid, "Z_MIN");
    float_type zmax  =  load_namedfloat(fid, "Z_MAX");
    char znettype[50];
    load_namedstringn(fid, "ZNET_TYPE", znettype, 50);
    ZNET  =  (float_type*)malloc(sizeof(float_type)*N_Z);
    create_net(zmin, zmax, N_Z, znettype, ZNET);

    ZSTEP = ZNET[1]-ZNET[0];
    CURRENT_Z = ZNET[0];		

}

void  create_mystartcondition(FILE* fid)
{
    #ifndef _SILENCE
	 printf("\n[%d]: Creating starting condition...", PROCESS_RANK);
    #endif
    char pulseshape[200] = ""; 
    load_namedstringn(fid, "PULSE_SHAPE", pulseshape, 200);
   
    if (ISMASTER)
    {
      float_type tau_fwhm = load_namedfloat(fid,   "DURATION_FWHM"); 
      float_type d_fwhm   = load_namedfloat(fid,   "DIAMETER_FWHM");
      float_type E        = load_namedfloat(fid,   "ENERGY");
      float_type f        = load_namedfloat(fid,   "FOCUSING_DISTANCE");
      float_type noiselevel = load_namedfloat(fid, "NOISE_LEVEL");
      printf("\npulse shape = %s", pulseshape);
       
      printf("\nBeam diameter FWHM  = %e", (double)d_fwhm);
      printf("\nPulse energy        = %e", (double)E);
      printf("\nFocusing distance   = %e", (double)f);
      printf("\nnoiselevel          = %e", (double)noiselevel);
      fflush(stdout);
    }
    
    if      (strncmp(SHAPE_GG, pulseshape,strlen(SHAPE_GG)) == 0) 									create_gg(fid);
}



void load_nonlindata(FILE* fid)
{
	NONLIN_REFRINDEX 	= load_namedfloat(fid, "NONLIN_REFRINDEX");

	float_type d = load_namedfloat(fid, "QUADRATIC_NONLINEARITY");
	d *= sqrt(2.0*OMEGA0/real(WAVENUMBER0[0])*VACUUM_PERMEABILITY)/2.0;

	QUADRATIC_NONLINEARITY_EEO = d*sin(2*THETA_OA)*cos(2*PHI_OA);
	QUADRATIC_NONLINEARITY_EOO = d*sin(THETA_OA)  *sin(2*PHI_OA);

	RAMAN_FRACTION   	= load_namedfloat(fid, "RAMAN_FRACTION");
	TAU_RAMAN 	 	= load_namedfloat(fid, "TAU_RAMAN");
	OMEGA_RAMAN      	= load_namedfloat(fid, "OMEGA_RAMAN");

	NEUTRAL_DENSITY    	    = load_namedfloat(fid, "NEUTRAL_DENSITY");
	RECOMBINATION_TAU  	    = load_namedfloat(fid, "RECOMBINATION_TAU");
	COLLISION_TAU        	= load_namedfloat(fid, "COLLISION_TAU");
	IONIZATION_POTENTIAL 	= load_namedfloat(fid, "IONIZATION_POTENTIAL");

#ifdef MULTIPHOTON_IONIZATION
	BETA_MPI		= load_namedfloat(fid, "MPI_CROSSSECTION");
#endif
#ifdef TUNNEL_IONIZATION
	TUNNELING_FIELD = load_namedfloat(fid, "TUNNELING_FIELD");
#endif


	IONIZATION_POTENTIAL *= ELECTRON_CHARGE;
	AMBIENT_CARRIER_DENSITY = load_namedfloat(fid, "AMBIENT_CARRIERS",true,0)*NEUTRAL_DENSITY;
}

float_type load_namedfloat(FILE* fid,const char* name, bool defvalue_present, float_type defvalue)
{
    //first, look for the parameter in command line.   
    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0 ) return (float_type)atof(_ARGV[i+1]);
    
    //then, in the input file.
	fseek(fid, 0, SEEK_SET);
	char namebuf[500];
	char valbuf[300];
	char buf[800];
	while (!feof(fid))
	{
		fgets(buf, 800, fid);
		sscanf(buf, "%s %s",namebuf, valbuf);
		if (strcmp(namebuf, name) == 0) return (float_type)atof(valbuf);
	}
	if (!defvalue_present)
	{
		printf("\nload_namedparameter: Paramneter %s not found.", name);
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
	char buf[800];
    while (!feof(fid))
    {
	 fgets(buf, 800, fid);
     sscanf(buf, "%s %s",namebuf, valbuf);
     if (strcmp(namebuf, namenum) == 0) return atof(valbuf);
    }

    for(int i=2; i<_ARGC-1; i++) if (strcmp(_ARGV[i],name)==0 ) return atof(_ARGV[i+1]);
    //then, in the input file.
    fseek(fid, 0, SEEK_SET);
    while (!feof(fid))
    {
	 fgets(buf, 800, fid);
     sscanf(buf, "%s %s",namebuf, valbuf);
     if (strcmp(namebuf, name) == 0) return atof(valbuf);
    }
    if (!defvalue_present)
    {
     printf("\nload_namedparameter: Paramneter %s not found.", name);
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
	char buf[800];
	while (!feof(fid))
	{
		fgets(buf, 800, fid);
		sscanf(buf, "%s %s",namebuf, valbuf);
		if (strncmp(namebuf, name, 500) == 0) return atoi(valbuf);
	}
	if (!defvalue_present)
	{
		printf("\nload_namedint: Paramneter %s not found.", name);
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
	char buf[800];
	while (!feof(fid))
	{
		fgets(buf, 800, fid);
		sscanf(buf, "%s %s",namebuf, valbuf);
		if (strncmp(namebuf, name, 500) == 0) return atoi(valbuf);
	}
	if (!defvalue_present)
	{
		printf("\nload_namedint: Paramneter %s not found.", name);
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
	char buf[800];
	while (!feof(fid))
	{
		fgets(buf, 800, fid);
		sscanf(buf, "%s",namebuf);
		if (strncmp(namebuf, name, 500) == 0) 
		{
			buf[strlen(buf)-1]=0;
			char* b = buf + strlen(namebuf); while ((*b) == ' ') b++;
			strncpy(output, b, N); 
		    return;
		}
	}
	if (!defvalue_present)
	{
		printf("\nload_namedstringn: Paramneter %s not found.", name);
		throw "load_namedstringn: parameter not found";
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
	char buf[800];
	while (!feof(fid))
	{
		fgets(buf, 800, fid);
		sscanf(buf, "%s",namebuf);
		if (strncmp(namebuf, name, 500) == 0) 
		{
			buf[strlen(buf)-1]=0;
			char* b = buf + strlen(namebuf); 
	        while ((*b)==' ') b++;
			strncpy(output, b, N); return;
		}
	}
	if (!defvalue_present)
	{
		printf("\nload_nameddouble: Paramneter %s not found.", name);
		throw "load_nameddouble: parameter not found";
	}
	else strncpy(output, defvalue, N); 
}



void print_variables()
{
    printf(  "\nTIME_START  =%e", (double)TIME_START);

#ifdef _FULL_DUMPS
    printf(  "\nDUMP_ID     =%d", DUMP_ID);
    printf(  "\nDUMP_PREFIX =%s", DUMP_PREFIX);
#endif

    printf("\n\nOMEGA0=%e", (double)OMEGA0);
    
    printf("\n\nTMIN=%e, TMAX=%e, N_T=%d", (double)TMIN, (double)TMAX, N_T);

    printf(  "\nN_X=%d, XMIN=%e, XMAX=%e",N_X, (double)XMIN, (double)XMAX);
    printf(  "\nN_Y=%d, YMIN=%e, YMAX=%e",N_Y, (double)YMIN, (double)YMAX); 
    
    /*printf(  "\nN_Z=%d, ZNET=",N_Z);         for (int i=0;i<N_Z;i++) printf(" %e", ZNET[i]);
    
    printf("\n \nOMEGA                      WAVENUMBER             ");
    for(int i=0;i<N_T;i++)
    {
     printf(  "\n%10e       %10e+i%10e", OMEGA[i], real(WAVENUMBER[i]), imag(WAVENUMBER[i]));
    }*/
	printf("\nN_Z=%d, ZNET[0]=%e, ZNET[%d]=%e", N_Z, ZNET[0], N_Z-1, ZNET[N_Z-1]);
	printf("\nOMEGA0 = %e, WAVENUMBER0o = %e + %ei", (double)OMEGA0, (double)real(WAVENUMBER0[0]), (double)imag(WAVENUMBER0[0]));
	printf("\nGROUP_VELOCITY=%e", (double)GROUP_VELOCITY);
	printf("\nQUADRATIC_NONLINEARITY_EEO=%e", (double)QUADRATIC_NONLINEARITY_EEO);
	printf("\nQUADRATIC_NONLINEARITY_EOO=%e", (double)QUADRATIC_NONLINEARITY_EOO);
    
	printf("\nTHETA_OA=%e", (double)THETA_OA);
	printf("\nPHI_OA=%e", (double)PHI_OA);
    
	
	printf("\n\nNONLIN_REFRINDEX = %e", (double)NONLIN_REFRINDEX);
    printf(  "\nRAMAN_FRACTION   = %e", (double)RAMAN_FRACTION);
    printf(  "\nTAU_RAMAN        = %e", (double)TAU_RAMAN);
    printf(  "\nOMEGA_RAMAN      = %e", (double)OMEGA_RAMAN);
    printf("\n\nNEUTRAL_DENSITY        = %e", (double)NEUTRAL_DENSITY);
    printf(  "\nAVALANCHE_CROSSSECTION = %e", (double)AVALANCHE_CROSSSECTION);
    printf(  "\nRECOMBINATION_TAU      = %e", (double)RECOMBINATION_TAU);
    printf(  "\nCOLLISION_TAU          = %e", (double)COLLISION_TAU);
    printf(  "\nIONIZATION_POTENTIAL   = %e", (double)IONIZATION_POTENTIAL);


#ifdef MULTIPHOTON_IONIZATION
    printf(  "\nBETA_MPI               = %e", (double)BETA_MPI);
    printf(  "\nK_MPI                  = %d", (double)K_MPI);
#endif
#ifdef TUNNEL_IONIZATION
	printf(  "\nTUNNELING_FIELD   = %e", (double)TUNNELING_FIELD);
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
/*	if (omega0 > omega_hw) 
	{
		for (int i=0;     i<Nt/2; i++) omega[i] = omega0 + wstep*i;
		for (int i=Nt/2;  i<Nt;   i++) omega[i] = omega0 - wstep*(Nt-i);
	}
	else
	{
		for (int i=0;     i<Nt; i++) omega[i] = wstep*(i+1);
	}
*/	
	for (int i=0;     i<Nt/2; i++) omega[i] = omega0 + wstep*i;
	for (int i=Nt/2;  i<Nt;   i++) omega[i] = omega0 - wstep*(Nt-i);

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
		if (!fid) {perror(""); throw "Unable to open file with refractive index information!";}
		fread(&filetype, sizeof(int), 1, fid);
	}
	f_complex* refr_index = (f_complex*)malloc_ch(2*sizeof(f_complex)*N);
	if (filetype < FILETYPE_SELLMEIER_LAMBDA && ISMASTER) printf("\n Warning! Refractive index file is a scalar one, assuming ne=no.");
	switch (filetype)
	{
		case FILETYPE_NOFILE           : for (int i=0; i<N; i++) refr_index[i] = 1.0;                         break;
		case FILETYPE_SELLMEIER_LAMBDA : load_refrindex_sellmeier_lambda(fid, N, omega, refr_index); break;
		case FILETYPE_SELLMEIER_OMEGA  : load_refrindex_sellmeier_omega (fid, N, omega, refr_index); break;
		case FILETYPE_RAW              : load_refrindex_raw             (fid, N, omega, refr_index); break;
		case FILETYPE_RAW_VEC		   : load_refrindex_raw_vec			(fid, N, omega, refr_index); break; 
		default: throw "Unknown type of file with refractive index information!";
	}

    //propagation direction correction:
	for (int i=0; i<N; i++) 
	{
		f_complex no = refr_index[2*i];
		f_complex ne = refr_index[2*i+1]; 
		refr_index[2*i+1] = ((float_type)1.0)/sqrt((float_type)sin(THETA_OA)*(float_type)sin(THETA_OA)/ne/ne + (float_type)cos(THETA_OA)*(float_type)cos(THETA_OA)/no/no);
	}

	f_complex j = f_complex(0,1);

	for (int i=0; i<N; i++) wavenum[2*i]    = conj(refr_index[2*i])  *f_complex(omega[i]/LIGHT_VELOCITY);
    for (int i=0; i<N; i++) wavenum[2*i+1]  = conj(refr_index[2*i+1])*f_complex(omega[i]/LIGHT_VELOCITY);

	free(refr_index);
	if (filetype != FILETYPE_NOFILE) fclose(fid);
}



void load_refrindex_sellmeier_omega (FILE* fid, int N, float_type* cfreq, f_complex* refrindex)
{
	int koefN = 0;
	double* sB = NULL;
	double* somega = NULL;
	fread(&koefN, sizeof(int), 1, fid);

	sB      = (double*)malloc_ch(koefN*sizeof(double));
	somega  = (double*)malloc_ch(koefN*sizeof(double));
	fread(somega,  sizeof(double), koefN, fid);
	fread(sB,      sizeof(double), koefN, fid);
		
	for (int i=0; i<N; i++)
	{
		f_complex n2 = 1.0;
		for (int j=0; j<koefN; j++) n2 += sB[j]*(somega[j]*somega[j])/(somega[j]*somega[j] - cfreq[i]*cfreq[i]);
		refrindex[2*i] = sqrt(n2); refrindex[2*i+1] = sqrt(n2);
	}
	delete sB;
	delete somega;
}

void load_refrindex_sellmeier_lambda(FILE* fid, int N, float_type* cfreq, f_complex* refrindex)
{
	int koefN = 0;
	double* sB	= NULL;
	double* slambda = NULL;
	double* somega  = NULL;
	fread(&koefN, sizeof(int), 1, fid);

	sB      = (double*)malloc_ch(koefN*sizeof(double));
	somega  = (double*)malloc_ch(koefN*sizeof(double));
	slambda = (double*)malloc_ch(koefN*sizeof(double));

	fread(slambda,  sizeof(double), koefN, fid);
	fread(sB,      sizeof(double), koefN, fid);
		
	for (int j=0; j<koefN; j++) somega[j] = 2*M_PI*LIGHT_VELOCITY/slambda[j];

	for (int i=0; i<N; i++)
	{
		f_complex n2 = 1.0;
		for (int j=0; j<koefN; j++) n2 += sB[j]*(somega[j]*somega[j])/(somega[j]*somega[j] - cfreq[i]*cfreq[i]);
		refrindex[2*i] = sqrt(n2);	refrindex[2*i+1] = sqrt(n2);
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
					float_type omega  = (float_type)(pomega[j]);
					float_type omega_ = (float_type)(pomega[j-1]);
                                        f_complex nj  = f_complex((float_type)pn[2*j],  (float_type)pn[2*j+1]);
					f_complex nj_ = f_complex((float_type)pn[2*j-2],(float_type)pn[2*j-1]);
					refrindex[2*i]   = nj_+(nj-nj_)*f_complex((fabs(cfreq[i])-omega_)/(omega-omega_));
					refrindex[2*i+1] = refrindex[2*i];
					break;
				}
			if (j==(pointsN-1))
			{
				refrindex[2*i]=1;
				refrindex[2*i+1]=1;
			}
		}
	}

    free(pomega);
	free(pn);
}

void load_refrindex_raw_vec(FILE* fid, int N, float_type* cfreq, f_complex* refrindex)
{
#ifdef _DEBUG
	printf("\n void load_refrindex_raw(FILE*, int, float_type*, float_type*)");
	printf("\n N = %d", N);
#endif
	int pointsN = 0;
	double*      pomega = NULL;
	double*      pn     = NULL;
	fread(&pointsN,sizeof(int),1, fid);

	pomega = (double*)malloc_ch(  sizeof(double)*pointsN);
	pn     = (double*)malloc_ch(4*sizeof(double)*pointsN);
	
	fread(pomega,sizeof(double),   pointsN, fid);
	fread(pn    ,sizeof(double), 4*pointsN, fid);


	OMEGA_MAX = pomega[pointsN-1]; OMEGA_MIN = pomega[0];
	for (int j=0; j<pointsN; j++) {OMEGA_MAX = max(OMEGA_MAX,pomega[j]); OMEGA_MIN = min(OMEGA_MIN, pomega[j]);}
	for (int i=0; i<N; i++)
	{
		for (int j=1; j<pointsN; j++)
		{
			if ((pomega[j-1]-cfreq[i])*(pomega[j]-cfreq[i]) <= 0) 
				{
					float_type omega  = (float_type)(pomega[j]);
					float_type omega_ = (float_type)(pomega[j-1]);

                    f_complex njo  = f_complex((float_type)pn[4*j]  ,  (float_type)pn[4*j+1]);
					f_complex nje  = f_complex((float_type)pn[4*j+2],  (float_type)pn[4*j+3]);
					f_complex njo_ = f_complex((float_type)pn[4*j-4],(float_type)pn[4*j-3]);
					f_complex nje_ = f_complex((float_type)pn[4*j-2],(float_type)pn[4*j-1]);

					refrindex[2*i]     = njo_+(njo-njo_)*f_complex((fabs(cfreq[i])-omega_)/(omega-omega_));
					refrindex[2*i+1]   = nje_+(nje-nje_)*f_complex((fabs(cfreq[i])-omega_)/(omega-omega_));

					break;
				}
			if (j==(pointsN-1))
			{
				refrindex[2*i]=1;
				refrindex[2*i+1]=1;
			}
		}
	}

    free(pomega);
	free(pn);
}

float_type calculate_groupvelocity(int Nt, float_type* omega, f_complex* wavenum, float_type omega0)
{
	float_type wstep = fabs(omega[1]-omega[0]);

	for (int i=0; i<Nt-1; i++) if (fabs(omega0-omega[i])<wstep) return (omega[i]-omega[i+1])/real(wavenum[2*i]-wavenum[2*i+2]);
	throw "calculate_groupvelocity: invalid input: omega0 is not inside the omega net";	

}


void init_zstep_kerr()
{
 float_type maxI = 0;
 for (int ny=0; ny<MY_NY; ny++)
 for (int nx=0; nx<N_X;   nx++)
 {
  for (int nt=0; nt<2*N_T; nt++)
  {
   int ofs =nt+2*N_T*(nx+N_X*ny);
   maxI = max(maxI, abs2(BIGBUFFER1[ofs]));
  }
 }
 float_type maxI_ = maxI;
 MPI_Allreduce(&maxI_,  &maxI,  1, MPI_FLOAT_TYPE, MPI_MAX, MPI_COMM_WORLD);  

 float_type kNL = OMEGA0/LIGHT_VELOCITY*NONLIN_REFRINDEX*maxI;
 if (ZSTEP*kNL > MAX_TOLERANCE) ZSTEP = MAX_TOLERANCE/kNL;
 if (ISMASTER) {printf("\n Initializing ZSTEP according to Kerr nonlinearity. maxI=%e, kNL=%e, ZSTEP=%e", maxI, kNL, ZSTEP); fflush(stdout);}
}




