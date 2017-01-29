#include <ionization.h>

void calculate_plasmadensity_small(f_complex* field, f_complex* pro)
{
	pro[0] = AMBIENT_CARRIER_DENSITY*2*frand();
	float_type tstep = (TMAX-TMIN)/N_T;
	for (int nt=1; nt<N_T;   nt++)
	{
	    // Plasma density calculation in performed via Runge-Khutta method of second order (Heun method)
	   
        float_type reE0 = real(field[nt-1]), imE0 = imag(field[nt-1]);
        float_type reE1 = real(field[nt])  , imE1 = imag(field[nt]); 
        float_type ro0  = real(pro[nt-1]);
        
        float_type k1 = (tstep*IONRATE_DENOM)*plasma_source_function(reE0, imE0, ro0);	    
        float_type k2 = (tstep*IONRATE_DENOM)*plasma_source_function(reE1, imE1, ro0+k1);
	    pro[nt]   = ro0 + 0.5*(k1+k2);
	}

}

void calculate_plasmadensity_small_2float(f_complex* field, float_type* pro)
{
	pro[0] = AMBIENT_CARRIER_DENSITY*2*frand();
	float_type tstep = (TMAX-TMIN)/N_T;
	for (int nt=1; nt<N_T;   nt++)
	{
	    // Plasma density calculation in performed via Runge-Khutta method of second order (Heun method)
	   
        float_type reE0 = real(field[nt-1]), imE0 = imag(field[nt-1]);
        float_type reE1 = real(field[nt])  , imE1 = imag(field[nt]); 
        float_type ro0  = pro[nt-1];
        
        float_type k1 = (tstep*IONRATE_DENOM)*plasma_source_function(reE0, imE0, ro0);	    
        float_type k2 = (tstep*IONRATE_DENOM)*plasma_source_function(reE1, imE1, ro0+k1);
	    pro[nt]   = ro0 + 0.5*(k1+k2);
	}
}

void calculate_plasmadensity(f_complex* field, f_complex* pro)
{
    
    for (int ny=0; ny<MY_NY; ny++)
	for (int nx=0; nx<N_X; nx++) 
	{
		int ofs = N_T*(nx+N_X*ny); calculate_plasmadensity_small(field + ofs, pro+ofs);
	}
}

void calculate_plasmadensity_2float(f_complex* field, float_type* pro)
{
    
    for (int ny=0; ny<MY_NY; ny++)
	for (int nx=0; nx<N_X; nx++) 
	{
		int ofs = N_T*(nx+N_X*ny); calculate_plasmadensity_small_2float(field + ofs, pro+ofs);
	}
}

float_type photoionization_function(float_type reA, float_type imA)
{
        float_type I = abs2(reA, imA);
        if (I < exp(IONIZATION_MIN_I_LN)) return 0;
	float_type lnI = log(I);

        //if (isnan(lnI)) throw "float_type photoionization_function(...): Intensity is nan!";

#ifdef DIRECT_IONIZATION_RATE
        float_type lnW = photoionization_rate_ln(lnI)-IONRATE_DENOM_LN;
	return exp(lnW);
#endif	

        int n = (int)floor((lnI - IONIZATION_MIN_I_LN)/IONIZATION_I_LN_TOLERANCE);
        if (n > IONIZATION_N-2) throw "Attention!!!! Intensity is too big! Keep away, blast is possible."; 
        float_type lnIn  = IONIZATION_MIN_I_LN+n*IONIZATION_I_LN_TOLERANCE;

        float_type lnWn  = IONIZATION_RATE_LN[n];
        float_type lnWnp = IONIZATION_RATE_LN[n+1];
        float_type lnW = lnWn + (lnI-lnIn)*(lnWnp-lnWn)/IONIZATION_I_LN_TOLERANCE;

	return exp(lnW);
}


float_type keldysh_rate_ln(float_type lnI)
{
	const int sumN = 30;
	float_type E = exp(lnI/2)*sqrt(2*OMEGA0/real(WAVENUMBER0[0])*VACUUM_PERMEABILITY);
	float_type m = 0.64*ELECTRON_MASS;
	float_type kgamma = OMEGA0*sqrt(m/ELECTRON_CHARGE/ELECTRON_CHARGE*IONIZATION_POTENTIAL)/E;
	float_type Gamma = kgamma*kgamma/(1+kgamma*kgamma);
	float_type Xi    = 1/(1+kgamma*kgamma);

	float_type Kg, Eg; ellipke(Gamma, &Kg, &Eg);
	float_type Kx, Ex; ellipke(Xi   , &Kx, &Ex);

	float_type alpha = M_PI*(Kg-Eg)/Ex, beta = M_PI*M_PI/4/Kx/Ex;
	
	float_type x     = 2/M_PI*IONIZATION_POTENTIAL/PLANCK_CONSTANT_REDUCED/OMEGA0*Ex/sqrt(Gamma); 
	float_type nu    = floor(x+1) - x;
	
	float_type Q = 0; for (int i=0;i<sumN;i++) Q+=exp(-i*alpha)*dawson(sqrt(beta*(i+2*nu)));
	Q=Q*sqrt(M_PI/2/Kx);
	
//	float_type W = 2*OMEGA0/9/M_PI*pow(OMEGA0*m/PLANK_CONSTANT_REDUCED/sqrt(Gamma), 1.5)*Q*exp(-alpha*floor(x+1));

	float_type lnW = log(2*OMEGA0/9/M_PI) + 1.5*log(OMEGA0*m/PLANCK_CONSTANT_REDUCED)-0.75*log(Gamma)+log(Q)-alpha*floor(x+1);
	return lnW;
}


void calc_ionization_rate()
{
	#ifndef _SILENCE
	printf("\n[%d]: Calculating inonization profile...",PROCESS_RANK); fflush(stdout);
	#endif
	
	for (int i=0; i<IONIZATION_N; i++)
	{
		float_type lnI    = IONIZATION_MIN_I_LN + IONIZATION_I_LN_TOLERANCE*i;
		IONIZATION_RATE_LN[i] = photoionization_rate_ln(lnI) - IONRATE_DENOM_LN;
	}	
	#ifndef _SILENCE
	printf("Done.");
	#endif
}





















