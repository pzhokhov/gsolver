#include "ionization.h"
#include <gsl/gsl_sf_gamma.h>


void load_ionization_profile(char* filename);

inline float_type arcsh(float_type x) {	return log(x+sqrt(x*x+1));}

void calculate_plasmadensity_losses_small(f_complex* field, float_type* pro, int stride, f_complex* loss, f_complex* n2factor)
{
	float_type tstep = TSTEP*IONRATE_DENOM;
	if (loss) loss[0]=0;
	if (n2factor) n2factor[0]=1.0; 
#ifndef MULTI_LEVEL_IONIZATION
	float_type ro1 = AMBIENT_CARRIER_DENSITY*2.0*frand();
	float_type ro2 = 0;
	float_type k1, k2;
	if (pro) pro[0] = ro1;
#else
	float_type* ro = ION_DENSITIES_BUFFER; float_type* W = IONIZATION_RATES_BUFFER;
	int N = IONIZATION_LEVEL_N; 
	for (int n=0; n<N; n++) ro[n]=0;
	if (pro) pro[0] = 0; 
#endif
	float_type reE0=0, imE0=0, reE1=0, imE1=0, I0=0, I1=0;

	for (int nt=1; nt<N_T;   nt++)
	{
	    // Plasma density calculation in performed via Runge-Khutta method of second order (Heun method)
	   	
	f_complex j = f_complex(0,1); 
	f_complex pfactor = exp(-j*OMEGA0*(TMIN + nt*(TMAX-TMIN)/N_T));
        float_type reE1 = real(field[nt]*pfactor), imE1 = imag(field[nt]*pfactor);
        float_type I1 = abs2(reE1, imE1);
	float_type tl=0;

#ifndef MULTI_LEVEL_IONIZATION
		float_type W01 = photoionization_function(reE0, imE0, ro1);    k1=(tstep)*(W01+avalanche_ionization_function(reE0, imE0, ro1)     - recombination_function(ro1));
		float_type W11 = photoionization_function(reE1, imE1, ro1+k1); k2=(tstep)*(W11+avalanche_ionization_function(reE1, imE1, ro1+k1)  - recombination_function(ro1+k1));

		ro1 += 0.5*(k1+k2);  if (ro1 > NEUTRAL_DENSITY) ro1 = NEUTRAL_DENSITY;	
		if (pro) pro[nt*stride] = ro1; 
		if (loss) loss[nt-1] = (I0>IONIZATION_MIN_I_LN)?(-(float_type)0.5*(W01)*(IONIZATION_POTENTIAL+PONDEROMOTIVE_COEFFICIENT*I0)/I0*field[nt-1]):(0);
		if (n2factor) n2factor[nt]=1.0;
#else
	 if (pro) pro[nt*stride]=pro[(nt-1)*stride];
	 if (loss) tl=0;
	 if (n2factor) n2factor[nt]=n2factor[nt-1];
 	 if (I0 > IONIZATION_MIN_I/INTENSITY_DENOM)
	 { 
	   photoionization_functionsN(reE0, imE0,     W);
	   photoionization_functionsN(reE1, imE1,   W+N);
	   float_type k1=0, k2=0;
#ifdef NO_UPPER_LEVEL_DEPLETION
	   if (IONIZATION_LEVEL_N == 1)
	   {k1 = tstep*W[0]*(NEUTRAL_DENSITY); k2 = tstep*W[N]*(NEUTRAL_DENSITY);}
	   else 
#endif
	   {k1 = tstep*W[0]*(NEUTRAL_DENSITY-ro[0]); k2 = tstep*W[N]*(NEUTRAL_DENSITY-ro[0]-k1);}

	   if (k1 < 0) k1=0; if (k2 < 0) k2 = 0;
	   ro[0] += 0.5*(k1+k2); 
#ifdef NO_UPPER_LEVEL_DEPLETION
	   if (IONIZATION_LEVEL_N > 1)
#endif
	   if (ro[0] > NEUTRAL_DENSITY) ro[0]=NEUTRAL_DENSITY; 
	   if (pro) pro[nt*stride] = ro[0]; 
           if (loss) tl = -(float_type)0.5*k1/tstep*(IONIZATION_POTENTIALS[0]+PONDEROMOTIVE_COEFFICIENT*I0)/I0;
	   if (n2factor) n2factor[nt] = 1.0-ro[0]/NEUTRAL_DENSITY*N2_IONFACTOR[0];

	   for (int n=1; n<IONIZATION_LEVEL_N; n++) 
	   {
#ifdef NO_UPPER_LEVEL_DEPLETION
	    if (n == IONIZATION_LEVEL_N-1) { k1 = tstep*W[n]*(ro[n-1]); k2 = tstep*W[N+n]*(ro[n-1]); }
            else
#endif
	    {  k1 = tstep*W[n]*(ro[n-1]-ro[n]); k2 = tstep*W[N+n]*(ro[n-1]-ro[n]-k1); }

	    if (k1 < 0) k1=0; if (k2 < 0) k2 = 0;
	    ro[n] = ro[n] + 0.5*(k1+k2); 
#ifdef NO_UPPER_LEVEL_DEPLETION
	    if (n < (IONIZATION_LEVEL_N-1)) 
#endif
	    if (ro[n]>ro[n-1]) ro[n]=ro[n-1];
	    if (pro) pro[nt*stride] += ro[n];
            if (loss) tl -= (float_type)0.5*k1/tstep*(IONIZATION_POTENTIALS[n]+PONDEROMOTIVE_COEFFICIENT*I0)/I0;
	    if (n2factor) n2factor[nt] -= ro[n]/NEUTRAL_DENSITY*N2_IONFACTOR[n]; 
	   } 
	 }
	 if (loss) loss[nt-1] = tl*field[nt-1];
#endif	  
	 reE0 = reE1; imE0 = imE1; I0=I1;
	}
	if (loss) loss[N_T-1]=0;
}


float_type photoionization_function(float_type reA, float_type imA, float_type ro)
{
        float_type I = abs2(reA, imA);
        if (I < (IONIZATION_MIN_I/INTENSITY_DENOM)) return 0;
		if (ro > NEUTRAL_DENSITY) return 0;
		float_type lnI = log(I)+INTENSITY_DENOM_LN;
        if (isnan(lnI)) 
			throw "float_type photoionization_function(...): Intensity is nan!";
#ifdef MULTIPHOTON_IONIZATION
	return exp(K_MPI*lnI+BETA_MPI_LN-IONRATE_DENOM_LN)*(NEUTRAL_DENSITY-ro);
#endif

#ifdef DIRECT_IONIZATION_RATE
        float_type lnW = photoionization_rate_ln(lnI)-IONRATE_DENOM_LN;
#else	
        int n = (int)floor((lnI - IONIZATION_MIN_I_LN)/IONIZATION_I_LN_TOLERANCE);
        if (n > IONIZATION_N-2) 
		{
			printf("\n IONIZATION_N = %d, lnI = %g, n = %d", IONIZATION_N, lnI, n);
			throw "Attention!!!! Intensity is too big! Keep away, blast is possible."; 
		}
		if (n < 0) return 0; 
        float_type lnIn  = IONIZATION_MIN_I_LN+n*IONIZATION_I_LN_TOLERANCE;

        float_type lnWn  = IONIZATION_RATE_LN[n];
        float_type lnWnp = IONIZATION_RATE_LN[n+1];
        float_type lnW = lnWn + (lnI-lnIn)*(lnWnp-lnWn)/IONIZATION_I_LN_TOLERANCE;
#endif
#ifdef IONIZATION_GAS
 	float_type F = 1; 
  #ifdef PPT_IONIZATION
  #ifdef YUDIN_IVANOV_CORRECTION
        float_type theta = atan(imA/reA);
	float_type phase_v = OMEGA0/real(WAVENUMBER0);
        float_type E = exp(lnI/2)*sqrt(2*VACUUM_PERMEABILITY*phase_v*phase_v/GROUP_VELOCITY);
        float_type g = sqrt(2.0*IONIZATION_POTENTIAL*INTENSITY_DENOM/IONRATE_DENOM/ELECTRON_CHARGE/ELECTRON_CHARGE*ELECTRON_MASS)*OMEGA0/E;

	if (g < MAX_YI_GAMMA)
        {
         float_type wl = OMEGA0*PLANCK_CONSTANT_REDUCED/HARTREE_ENERGY; 
         float_type T  = (E*E/ATOMIC_FIELD/ATOMIC_FIELD)/wl/wl/wl;
         F = exp(-T*(YI_Phi(theta, g)-YI_Phi(0,g))); 
        }
  #endif
  #endif
	return exp(lnW+log(NEUTRAL_DENSITY-ro))*F;
#else
	return exp(lnW);
#endif
}


#ifdef MULTI_LEVEL_IONIZATION

void photoionization_functionsN(float_type reA, float_type imA, float_type* W)
{

        float_type I = abs2(reA, imA);
        if (I < (IONIZATION_MIN_I/INTENSITY_DENOM)) {for (int i=0; i<IONIZATION_LEVEL_N; i++) W[i]=0; return;}

	float_type lnI = log(I)+INTENSITY_DENOM_LN;
        if (isnan(lnI)) 
			throw "float_type photoionization_functionsN(...): Intensity is nan!";
#ifdef DIRECT_IONIZATION_RATE
        W[0] = exp(photoionization_rateN_ln(lnI, 0)-IONRATE_DENOM_LN)*(NEUTRAL_DENSITY-ro[0]); 
        for (int i=1; i<IONIZATION_LEVEL_N; i++) W[i] = exp(photoionization_rateN_ln(lnI, i)-IONRATE_DENOM_LN)*(ro[i-1]-ro[i]); 
#else
         int n = (int)floor((lnI - IONIZATION_MIN_I_LN)/IONIZATION_I_LN_TOLERANCE);
        if (n > IONIZATION_N-2)
	 {
	     printf("\n IONIZATION_N = %d, lnI = %g, n = %d", IONIZATION_N, lnI, n);
             throw "photoionization_functionsN(): Attention!!!! Intensity is too big! Keep away, blast is possible."; 
	 }
 	if (n < 0) throw "photoionization_functionsN(): n<0";

#ifdef YUDIN_IVANOV_CORRECTION
         float_type theta = atan(imA/reA);  
         float_type wl = OMEGA0*PLANCK_CONSTANT_REDUCED/HARTREE_ENERGY; 
	 float_type phase_v = OMEGA0/real(WAVENUMBER0);
         float_type E = exp(lnI/2)*sqrt(2*VACUUM_PERMEABILITY*phase_v*phase_v/GROUP_VELOCITY);
         float_type T  = (E*E/ATOMIC_FIELD/ATOMIC_FIELD)/wl/wl/wl;
#endif
 	
	for (int i=0; i<IONIZATION_LEVEL_N; i++)
	{
         float_type lnIn  = IONIZATION_MIN_I_LN+n*IONIZATION_I_LN_TOLERANCE;
 	 float_type lnWn  = IONIZATION_RATE_LN[IONIZATION_N*(i)+n];
         float_type lnWnp = IONIZATION_RATE_LN[IONIZATION_N*(i)+n+1];
	 float_type lnW = lnWn + (lnI-lnIn)*(lnWnp-lnWn)/IONIZATION_I_LN_TOLERANCE;
 	 float_type F = 1; 
  #ifdef PPT_IONIZATION
  #ifdef YUDIN_IVANOV_CORRECTION
	float_type g = sqrt(2.0*IONIZATION_POTENTIALS[i]*INTENSITY_DENOM/IONRATE_DENOM/ELECTRON_CHARGE/ELECTRON_CHARGE*ELECTRON_MASS)*OMEGA0/E;

	if (g < MAX_YI_GAMMA) F = exp(-T*(YI_Phi(theta, g)-YI_Phi(0,g)));         
  #endif
  #endif
	 W[i] = exp(lnW)*F;
	}
#endif
}
#endif



float_type keldysh_rate_ln(float_type lnI)
{
	double k0 = real(WAVENUMBER0);
	double U = IONIZATION_POTENTIAL; 
	double phase_v = OMEGA0/k0;
	double E = exp(lnI/2)*sqrt(2*phase_v/GROUP_VELOCITY*phase_v*VACUUM_PERMEABILITY);
	double m = dELECTRON_MASS*load_namedfloat(SC_FID, "REDUCED_MASS", true, 1); 

	double kgamma = OMEGA0*sqrt(m/ELECTRON_CHARGE/ELECTRON_CHARGE*U)/E;
	double Gamma = kgamma*kgamma/(1+kgamma*kgamma);
	double Xi    = 1/(1+kgamma*kgamma);

	double Kg, Eg; ellipke(Gamma, &Kg, &Eg);
	double Kx, Ex; ellipke(Xi   , &Kx, &Ex);

	double alpha = M_PI*(Kg-Eg)/Ex, beta = M_PI*M_PI/4.0/Kx/Ex;
	
	double x     = 2/M_PI*U/PLANCK_CONSTANT_REDUCED/OMEGA0*Ex/sqrt(Gamma); 
	double nu    = floor(x+1) - x;

	int sumN = (int)ceil(10/alpha); sumN = min(sumN, 1000);
	double Q = 0; 
	for (int i=0;i<sumN;i++) Q+=exp(-i*alpha)*dawson(sqrt(beta*(i+2*nu)));
	Q=Q*sqrt(M_PI/2/Kx);
	
//	double W = 2*OMEGA0/9/M_PI*pow(OMEGA0*m/PLANK_CONSTANT_REDUCED/sqrt(Gamma), 1.5)*Q*exp(-alpha*floor(x+1));

	double lnW = log(2*OMEGA0/9/M_PI) + 1.5*log(OMEGA0*m/PLANCK_CONSTANT_REDUCED)-0.75*log(Gamma)+log(Q)-alpha*floor(x+1);
	return (float_type)lnW;
}


float_type PPT_rate_ln(float_type lnI, float_type Ui, float_type Z, int l)
{
#ifdef DIRECT_IONIZATION_RATE
	Ui /= IONRATE_DENOM/INTENSITY_DENOM; 	
#endif
	double k0 = real(WAVENUMBER0);
	double phase_v = OMEGA0/k0;
	double E = exp(lnI/2)*sqrt(2*dVACUUM_PERMEABILITY*phase_v*phase_v/GROUP_VELOCITY);
	double m = dELECTRON_MASS;

	double g = OMEGA0*sqrt(2*m/dELECTRON_CHARGE/dELECTRON_CHARGE*Ui)/E; 
    double Ueff = Ui*(1+1.0/2.0/g/g);
    double Uh = dHARTREE_ENERGY;
    
    double nu = Ueff/dPLANCK_CONSTANT_REDUCED/OMEGA0;
    int nu0 = (int)floor(1+nu);
    double alpha = 2*(arcsh(g)-g/sqrt(1+g*g));
    double beta = 2*g/sqrt(1+g*g);
    double E0 = sqrt(2*Ui/dELECTRON_CHARGE*Ui/dELECTRON_CHARGE*Ui/dPLANCK_CONSTANT_REDUCED*m/dPLANCK_CONSTANT_REDUCED);
    float neff = Z/sqrt(2*Ui/Uh); 
  
    double Q = 0; 
	if (g > 0.1)
	{
	 int sumN = (int)ceil(10/alpha); 
     for (int k=0; k<sumN; k++) 
	 {
		double P = exp(-alpha*(k+nu0-nu))*dawson(sqrt(beta*(k+nu0-nu)));
		Q += P;
	 }
	}
	else Q = sqrt(3*M_PI)/4/g/g;
	
    double Cnl2 = pow(2.0, (2.0*neff))/(neff*gsl_sf_gamma(2*neff));
    
    
    double lnW = log(4*sqrt(2.0)/M_PI*Cnl2) + log((2*E0/E/sqrt(1+g*g)))*(2*neff-3.0/2.0);
		       lnW += (-2*nu*(arcsh(g)-g*sqrt(1+g*g)/(1+2*g*g)));
			   lnW += log((2*l+1)*Ui/dPLANCK_CONSTANT_REDUCED*g*g/(1+g*g)*Q);
 
#ifdef YUDIN_IVANOV_CORRECTION
    if (g < MAX_YI_GAMMA)
    {
     size_t Nphi = 1024;
     double F=0;
     double wl = OMEGA0*dPLANCK_CONSTANT_REDUCED/dHARTREE_ENERGY; 
     double T = (E*E/dATOMIC_FIELD/dATOMIC_FIELD)/wl/wl/wl;
    
     for (int i=0; i<Nphi; i++)
     {
      double theta = i*M_PI/2.0/Nphi;
      F+=exp(-T * (YI_Phi(theta,g)-YI_Phi(0, g)));
     }     
     F*=M_PI/2.0/Nphi;
     lnW -= log(F);
    } 
#endif
  
    return (float_type)lnW;
}

#ifdef YUDIN_IVANOV_CORRECTION

float_type YI_Phi(float_type theta, float_type g)
{
 float_type sth = fabs(sin(theta)); 
 float_type sth2 = sth*sth; 

 float_type a = 1.0 + g*g - sth2;
 float_type b = sqrt(a*a + 4*g*g*sth2); 
 float_type c = sqrt(pow(sqrt((a+b)/2.0)+g, 2.0) + pow(sqrt((b-a)/2.0)+sth, 2.0)); 
 
 return ((g*g + sth2 + 0.5)*log(c) - (3.0*sqrt((b-a)/2.0)/2.0)*sth - sqrt((a+b)/2.0)/2.0*g); 
}

#endif


void calc_ionization_rate()
{
	char ionprofile_filename[512];
	load_namedstringn(SC_FID, "IONIZATION_PROFILE", ionprofile_filename, 512, true, "NO");
	if (strncmp(ionprofile_filename, "NO", 512) == 0)
	{
	 #ifndef _SILENCE
	 if (ISMASTER) printf("\nCalculating ionization profile..."); fflush(stdout);
	 #endif
	

	 for (int i=0; i<IONIZATION_N-1; i++)
	 {
		float_type lnI    = IONIZATION_MIN_I_LN + IONIZATION_I_LN_TOLERANCE*i;
		IONIZATION_RATE_LN[i] = photoionization_rate_ln(lnI) - IONRATE_DENOM_LN;
#ifdef MULTI_LEVEL_IONIZATION
		for (int n=1; n<IONIZATION_LEVEL_N; n++) IONIZATION_RATE_LN[IONIZATION_N*n+i] = photoionization_rateN_ln(lnI, n)-IONRATE_DENOM_LN;
#endif
	 }	

	 #ifndef _SILENCE
	 if (ISMASTER) printf("Done.");
	 #endif

	}
	else
	{
		load_ionization_profile(ionprofile_filename);
	}

#ifdef DUMP_IONRATE
	 if (ISMASTER)
	 {
	  FILE* fid = fopen("ionrate.txt", "w");
	  for (int i=0; i<IONIZATION_N; i++)
	  {
		float_type lnI    = IONIZATION_MIN_I_LN + IONIZATION_I_LN_TOLERANCE*i;
		fprintf(fid, "%e %e\n", lnI, IONIZATION_RATE_LN[i]);
	  }
	  fclose(fid);
	 }
#endif
}


void load_ionization_profile(char* filename)
{
	FILE* fid = fopen(filename, "rt"); if (!fid) throw "Unable to open ionization profile!";
	char buf[512];

	double I_, W_; int N=0;
	double *I, *W;
	fgets(buf, 512, fid);
	int ncols = sscanf(buf, "%lf %lf", &I_, &W_);  
	if (ncols != 2) throw "Incorrect ionization profile file - it should contain two columns - first for field intensity, second for ionization rate";
	while (!feof(fid)) {fgets(buf, 512, fid); N++;} rewind(fid);

	if (N != IONIZATION_N) throw "Incorrect file length! Revise ionization profile or ionization parameters";

	I = (double*)malloc(sizeof(double)*N);
	W = (double*)malloc(sizeof(double)*N);

	for (int i=0; i<N; i++) 
		fscanf(fid, "%lf %lf", I+i, W+i);
		
	fclose(fid);
	
	for (int i=0; i<N; i++) IONIZATION_RATE_LN[i] = log(W[i])-IONRATE_DENOM;

	free(I);
	free(W);
}


















