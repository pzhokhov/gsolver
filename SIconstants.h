#ifndef _SI_CONSTANTS
#define _SI_CONSTANTS
#include "extmath.h"
#define LIGHT_VELOCITY				   ((float_type)299792458.0)
const float_type VACUUM_PERMITTIVITY   = 8.854187817e-12;
const float_type VACUUM_PERMEABILITY   = 1.2566370615e-6;

const float_type ELECTRON_CHARGE       = 1.60217733e-19;
const float_type ELECTRON_MASS	       = 9.1093897e-31;
const float_type PLANCK_CONSTANT       = 6.6260755e-34;
const float_type PLANCK_CONSTANT_REDUCED =1.0545727e-34;
	
const float_type AVOGADRO_NUMBER       = 6.0221367e23;
const float_type HARTREE_ENERGY        = 4.359748124303985e-018;

const float_type VACUUM_IMPEDANCE      = sqrt(VACUUM_PERMITTIVITY/VACUUM_PERMEABILITY);
const float_type BOHR_RADIUS            = 5.2917720859e-11;
const float_type ATOMIC_FIELD           = 5.1422085489e+11; 
const float_type ATOMIC_FIELD_          = ATOMIC_FIELD/sqrt(2.0*VACUUM_IMPEDANCE); 
const float_type ATOMIC_TIME            = 2.418884326505e-17;

const float_type BOLTZMANN_K 		= 1.380656631799432e-23;


const double dLIGHT_VELOCITY			= 299792458.0;
const double dVACUUM_PERMITTIVITY		= 8.854187817e-12;
const double dVACUUM_PERMEABILITY		= 1.2566370615e-6;

const double dELECTRON_CHARGE			= 1.60217733e-19;
const double dELECTRON_MASS				= 9.1093897e-31;
const double dPLANCK_CONSTANT			= 6.6260755e-34;
const double dPLANCK_CONSTANT_REDUCED	= 1.0545727e-34;
	
const double dAVOGADRO_NUMBER			= 6.0221367e23;
const double dHARTREE_ENERGY			= 4.359748124303985e-018;

const double dATOMIC_FIELD		= 5.142208548937811e+11;
const double dBOLTZMANN_K		= 1.380656631799432e-23;


#endif
