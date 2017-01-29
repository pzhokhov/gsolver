#ifndef _IONIZATION_HEADER
#define _IONIZATION_HEADER

#include "solver.h"

const float_type IONRATE_DENOM    = 1e20;
const float_type IONRATE_DENOM_LN = log(IONRATE_DENOM);


float_type keldysh_rate_ln(float_type lnI);
float_type plasma_source_function(float_type reA, float_type imA, float_type ro);
float_type absorbtion_function(float_type reA, float_type imA, float_type ro);

float_type photoionization_function(float_type reA, float_type imA);
float_type photoabsorbtion_function(float_type reA, float_type imA);
float_type avalanche_ionization_function(float_type reA, float_type imA, float_type ro);
float_type recombination_function(float_type ro);



inline float_type photoionization_rate_ln(float_type lnI)
{
	#ifdef MULTIPHOTON_IONIZATION
		return K_MPI*lnI + log(BETA_MPI) + log(NEUTRAL_DENSITY);
	#endif
	#ifdef TUNNEL_IONIZATION
		throw "photoionization_rate_ln() : Tunnel ionization is not supported yet!";
	#endif
	#ifdef KELDYSH_IONIZATION
		return keldysh_rate_ln(lnI);
	#endif
}


inline float_type plasma_source_function(float_type reA, float_type imA, float_type ro)
{
	return  photoionization_function(reA, imA) + avalanche_ionization_function(reA, imA, ro) - recombination_function(ro);
}



inline float_type photoabsorbtion_function(float_type reA, float_type imA)
{
	return (abs2(reA, imA)==0)?(0):(0.5*photoionization_function(reA, imA)/abs2(reA, imA)*(IONIZATION_POTENTIAL*IONRATE_DENOM));
	
}

inline float_type avalanche_ionization_function(float_type reA, float_type imA, float_type ro)
{
	double I = abs2(reA, imA); 
	return (ro>0)?(AVALANCHE_CROSSSECTION/(IONIZATION_POTENTIAL+PONDEROMOTIVE_COEFFICIENT*I)*I*(ro/IONRATE_DENOM)):0;
}
inline float_type recombination_function(float_type ro)
{
	return (ro>0)?(ro/IONRATE_DENOM/RECOMBINATION_TAU):0;
}


#endif

