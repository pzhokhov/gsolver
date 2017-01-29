#ifndef _IONIZATION_HEADER
#define _IONIZATION_HEADER

#include "solver.h"


float_type keldysh_rate_ln(float_type lnI);
float_type PPT_rate_ln(float_type lnI, float_type Ui, float_type Z, int l);
float_type plasma_source_function(float_type reA, float_type imA, float_type ro);
float_type absorbtion_function(float_type reA, float_type imA, float_type ro);

#ifdef YUDIN_IVANOV_CORRECTION
float_type YI_Phi(float_type theta, float_type g);
#define MAX_YI_GAMMA (20)
#endif

float_type photoionization_function(float_type reA, float_type imA, float_type ro);

#ifdef MULTI_LEVEL_IONIZATION
void photoionization_functionsN(float_type reA, float_type imA, float_type* W);
#endif 
float_type photoabsorbtion_function(float_type reA, float_type imA, float_type ro);
float_type avalanche_ionization_function(float_type reA, float_type imA, float_type ro);
float_type recombination_function(float_type ro);

inline float_type photoionization_rate_ln(float_type lnI)
{
	#ifdef MULTIPHOTON_IONIZATION
		return 0; //photoionization is calculated in a direct fashion
	#endif
	#ifdef TUNNEL_IONIZATION
		throw "photoionization_rate_ln() : Tunnel ionization is not supported yet!";
	#endif
	#ifdef KELDYSH_IONIZATION
		return keldysh_rate_ln(lnI);
	#endif
	#ifdef PPT_IONIZATION
		return PPT_rate_ln(lnI, IONIZATION_POTENTIAL, load_namedfloat(SC_FID, "PPT_Z",true,1), load_namedint(SC_FID, "PPT_L",true,0));
	#endif 
}

#ifdef MULTI_LEVEL_IONIZATION

inline float_type photoionization_rateN_ln(float_type lnI, int level)
{
	#ifdef PPT_IONIZATION
	    return PPT_rate_ln(lnI, IONIZATION_POTENTIALS[level], level+1, load_namedint(SC_FID, "PPT_L",true,0));
	#endif 
}

#endif




inline float_type plasma_source_function(float_type reA, float_type imA, float_type ro)
{
	return  photoionization_function(reA, imA,ro) + avalanche_ionization_function(reA, imA, ro) - recombination_function(ro);
}



inline float_type photoabsorbtion_function(float_type reA, float_type imA, float_type ro)
{
	float_type I = abs2(reA, imA);
	return (I < IONIZATION_MIN_I/IONRATE_DENOM)?(0):(0.5*photoionization_function(reA, imA, ro)/I*(IONIZATION_POTENTIAL+PONDEROMOTIVE_COEFFICIENT*I));
}


inline float_type photoabsorbtion_function2(float_type reA, float_type imA, float_type ro1, float_type ro2)
{
	float_type I = abs2(reA, imA);
	return (I < IONIZATION_MIN_I/IONRATE_DENOM)?(0):(0.5*photoionization_function(reA, imA, ro1)/I*(IONIZATION_POTENTIAL+PONDEROMOTIVE_COEFFICIENT*I));
}


inline float_type avalanche_ionization_function(float_type reA, float_type imA, float_type ro)
{

	if (ro < 0 || ro > NEUTRAL_DENSITY) return 0;
	float_type I = abs2(reA, imA); 
	return AVALANCHE_CROSSSECTION/(IONIZATION_POTENTIAL+PONDEROMOTIVE_COEFFICIENT*I)*I*ro*(1-ro/NEUTRAL_DENSITY);
}
inline float_type recombination_function(float_type ro)
{
	return (ro>0)?(ro/RECOMBINATION_TAU):0;
}


#endif

