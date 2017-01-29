#include "solver.h"
#include "pulseshapes.h"

#include <gsl/gsl_sf_gamma.h>

void create_gg(FILE* fid)
{
        f_complex j = f_complex(0,1);
        float_type tau_fwhm = load_namedfloat(fid, "DURATION_FWHM", true, -1); tau_fwhm = load_namedfloat(fid, "DURATION_FWHM_OMEGA", true, tau_fwhm*OMEGA0)/OMEGA0; 
	if (tau_fwhm < 0) throw "create_gg(): DURATION_FWHM or DURATION_FWHM_OMGEA should be specified!";
    	float_type d_fwhm   = load_namedfloat(fid, "DIAMETER_FWHM");
    	float_type noiselevel = load_namedfloat(fid, "NOISE_LEVEL", true, 0);
	float_type noiselevel_uni = load_namedfloat(fid, "NOISE_LEVEL_UNIFORM",true,0);
		float_type spatialnoiselevel = load_namedfloat(fid, "NOISE_LEVEL_SPATIAL", true, 0); 
 	float_type ce_phase  = load_namedfloat(fid, "CE_PHASE", true, 0); 
	
	float_type sg_T  = load_namedfloat(fid, "PULSE_SG_ORDER", true, 1); 
	float_type sg_r  = load_namedfloat(fid, "BEAM_SG_ORDER",  true, 1); 

	float_type phase_flip        = load_namedfloat(fid, "PHASE_FLIP", true, 0)*M_PI; 
	float_type phase_dip	     = load_namedfloat(fid, "PHASE_DIP", true, 0);
  	float_type phase_flip_l      = load_namedfloat(fid, "PHASE_FLIP_LAMBDA", true, 1); 
	float_type phase_flip_width  = load_namedfloat(fid, "PHASE_FLIP_WIDTH",  true, 0); 

	float_type phase_flip_maxw = 2*M_PI*LIGHT_VELOCITY/(phase_flip_l - phase_flip_width/2.0); 
	float_type phase_flip_minw = 2*M_PI*LIGHT_VELOCITY/(phase_flip_l + phase_flip_width/2.0); 


     	if (ISMASTER) printf("\nCreating super-gaussian-shape pulse with pulse order %g and beam order %g...",sg_T, sg_r );
	
 //       float_type A0 = sqrt(E/d_fwhm/d_fwhm/tau_fwhm*sqrt(log(16.0)/M_PI)*log(16.0)/M_PI);              

    	float_type E        = load_namedfloat(fid, "ENERGY", true, -1); 
	float_type E_P_factor = 2*pow(log(2), 1.0/2.0/sg_T)/gsl_sf_gamma(1.0/2.0/sg_T);
	float_type Pcr = 2*M_PI*LIGHT_VELOCITY/real(WAVENUMBER0)/OMEGA0/NONLIN_REFRINDEX*INTENSITY_DENOM; 
        float_type P0=load_namedfloat(fid, "POWER_CR", true, E/tau_fwhm*E_P_factor/Pcr)*Pcr; P0 = load_namedfloat(fid, "POWER", true, P0);
	E = P0*tau_fwhm/E_P_factor; 
	if (E<0) throw "create_gg(): ENERGY or POWER parameter should be specified!";

	float_type A0 = sqrt(P0/d_fwhm/d_fwhm*4*pow(log(2), 1.0/sg_r)/M_PI/gsl_sf_gamma(1.0/sg_r));

        float_type tau = tau_fwhm/2.0/pow(log(sqrt(2.0)), 1.0/2.0/sg_T);
        float_type ro  = d_fwhm  /2.0/pow(log(sqrt(2.0)), 1.0/2.0/sg_r);
        float_type zR  = real(WAVENUMBER0)*ro*ro/2;
            
	float_type C  = load_namedfloat(fid, "CHIRP_SI", true, 0);
	           C  = load_namedfloat(fid, "CHIRP",true,C*tau*tau);
		float_type Cw = load_namedfloat(fid, "CHIRP_SPECTRAL", true, 0);
	
		float_type f  = load_namedfloat(fid, "FOCUSING_DISTANCE", true, 0);
		f = load_namedfloat(fid, "FOCUSING_DISTANCE_DIFFRACTION", true, f/zR); 
        if (zR > 2*f) f = zR*zR  /2/f*(1-sqrt(1-4*f*f/zR/zR)); //inital wavefront curvature radius calculated from distance from linear focus
		f        = load_namedfloat(fid, "WAVEFRONT_RADIUS", true, f);	       
		f        = load_namedfloat(fid, "WAVEFRONT_RADIUS_DIFFRACTION", true, f/zR)*zR;      
	
	float_type zmin = load_namedfloat(fid, "Z_MIN_DIFFRACTION", true, ZNET[0]/zR)*zR; 
	float_type zmax = load_namedfloat(fid, "Z_MAX_DIFFRACTION", true, ZNET[N_Z-1]/zR)*zR; 
  	char znettype[50]; load_namedstringn(fid, "ZNET_TYPE", znettype, 50);
    	create_net(zmin, zmax, N_Z, znettype, ZNET);	
	            
        if (ISMASTER)
        {
         printf("\n   A0  = %e", A0);
	 printf("\n   P0  = %e", P0);
	 printf("\n   E   = %e", E);
         printf("\n   tau = %e", tau);
	 printf("\n   T_FWHM = %e", tau_fwhm); 
         printf("\n   ro  = %e", ro);
	 printf("\n   d_FWHM = %e", d_fwhm); 
	 printf("\n   C   = %e", C);
	 printf("\n   Cw  = %e", Cw);
	 printf("\n   f   = %e", f);
	 printf("\n   zR  = %e", zR); 
	 printf("\n   zmin = %e", zmin);
	 printf("\n   zmax = %e", zmax);
	 printf("\n   CEP = %g \n", ce_phase);
	
        }
                        
	for (int cy=0; cy<MY_NY; cy++)
	for (int cx=0; cx<MY_NX; cx++)
	{
	 int ofs0 = N_T*(cx+MY_NX*cy);                
     float_type R=0; 

#if     TRANSVERSE_DIMENSIONS == 1
#ifndef _UNIAXIAL_FINITE_DIFFERENCE
	 R = HT_PLAN->x_n(cx+MY_NXstart)*XMAX;
#else
	 R = XMAX*(cx+MY_NXstart)/N_X;
#endif
#elif   TRANSVERSE_DIMENSIONS == 2
	 float_type x = XMIN+cx*XSTEP, y = YMIN+(cy+MY_NYstart)*YSTEP;
	 R = sqrt(x*x+y*y);
#endif

	 f_complex Amax = A0*(float_type)(exp(-pow(R/ro, 2.0*sg_r)))*frand1c(spatialnoiselevel); //multiply this expression by exp(j*phase_function(x, y))
         for (int ct=0; ct<N_T; ct++) 
         {
	  float_type t = (TMIN + TSTEP*ct);
	  float_type T = t/tau; 
	  float_type carrier_phase = ce_phase - T*T*C + ((OMEGA0 == OMEGA[0])?(0):(-OMEGA0*(TMIN+ct*(TMIN-TMAX)/N_T)));
			
	  FIELD_REAL_SMALL[ct]   =   Amax*(float_type)exp(-(pow(T, 2.0*sg_T)))*(exp(j*carrier_phase)+noiselevel*(frand0()+j*frand0()))+A0*noiselevel_uni*(frand0()+j*frand0()); 
         }

         fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)FIELD_REAL_SMALL, (fftwt_complex*)(NL_SMALLBUFFER1)); 
	 for (int cw=0; cw<N_T; cw++)
	 {
#ifndef NO_SPACE_TIME_FOCUSING
	  float_type k = real(WAVENUMBER[cw]);
#else
	  float_type k = real(WAVENUMBER[0]);
#endif
	  float_type phase = (R*R*k/2.0/f) - Cw*(OMEGA[cw]-OMEGA0)*(OMEGA[cw]-OMEGA0)/2;
	
          if (phase_flip != 0 && phase_flip_minw < OMEGA[cw] && OMEGA[cw] < phase_flip_maxw) phase += phase_flip; 

	  FIELD[ofs0+cw]    = (NL_SMALLBUFFER1[cw]*exp(j*phase));
          if (phase_dip != 0 && phase_flip_minw < OMEGA[cw] && OMEGA[cw] < phase_flip_maxw) FIELD[ofs0+cw] *= (1.0-phase_dip); 
	 } 
        }
        if (ISMASTER) printf("Done.");
}
