#include "solver.h"
#include "pulseshapes.h"

void create_gg(FILE* fid)
{
        f_complex j = f_complex(0,1);
        float_type tau_fwhm = load_namedfloat(fid, "DURATION_FWHM"); 
    	float_type d_fwhm   = load_namedfloat(fid, "DIAMETER_FWHM");
    	float_type E        = load_namedfloat(fid, "ENERGY");
    	float_type f        = load_namedfloat(fid, "FOCUSING_DISTANCE");
    	float_type noiselevel = load_namedfloat(fid, "NOISE_LEVEL", true, 0);
		float_type noiselevel_uni = load_namedfloat(fid, "NOISE_LEVEL_UNIFORM",true,0);
		char polarization_str[128]; load_namedstringn(fid, "POLARIZATION",polarization_str,128); 
		f_complex pol[2]; sscanf(polarization_str, "%f %f %f %f",((float_type*)pol),((float_type*)pol)+1,((float_type*)pol)+2,((float_type*)pol)+3); 


		f_complex p0, p1;
		p0 = pol[0]/(float_type)sqrt((abs2(pol[0])+abs2(pol[1])));
		p1 = pol[1]/(float_type)sqrt((abs2(pol[0])+abs2(pol[1])));
		
     	if (ISMASTER) printf("\nCreating gaussian-shape pulse with polarization (%3f+i%3f, %3f+i%3f)",real(pol[0]),imag(pol[0]),real(pol[1]),imag(pol[1]));

        float_type A0 = sqrt(E/d_fwhm/d_fwhm/tau_fwhm*sqrt(log(16.0)/M_PI)*log(16.0)/M_PI);              
        float_type tau = tau_fwhm/sqrt(log(4.0));
        float_type ro  = d_fwhm/sqrt(log(4.0));
        //???
		float_type zR  = real(WAVENUMBER0[0])*ro*ro/2.0;
            
	
        if (zR > 2*f) f = zR*zR  /2/f*(1-sqrt(1-4*f*f/zR/zR));
			            
        if (ISMASTER)
        {
         printf("\n   A0  = %e", A0);
         printf("\n   tau = %e", tau);
         printf("\n   ro  = %e", ro);
        }
                        
	for (int cy=0; cy<MY_NY; cy++)
	for (int cx=0; cx<N_X; cx++)
	{
	 int ofs0 = 2*N_T*(cx+N_X*cy);    
	 float_type R; 
#ifdef _UNIAXIAL
	 R = FHT_PLAN->x_n(cx)*XMAX;
#else
	 float_type x = XMIN+cx*XSTEP, y = YMIN+(cy+MY_NYstart)*YSTEP;
	 R = sqrt(x*x+y*y);
#endif

	 float_type Amax = A0*(exp(-R*R/ro/ro)); 
     for (int ct=0; ct<N_T; ct++) 
     {
	  float_type t = (TMIN + TSTEP*ct);
	  float_type T = t/tau; 
	  float_type carrier_phase = (OMEGA0 == OMEGA[0])?(0):(-OMEGA0*(TMIN+ct*(TMIN-TMAX)/N_T));
			
	  FIELD_REAL_SMALL[2*ct]     =   p0*Amax*(float_type)exp(-T*T)*(exp(j*carrier_phase)+noiselevel*(frand0()+j*frand0()))+A0*noiselevel_uni*(frand0()+j*frand0()); 
	  FIELD_REAL_SMALL[2*ct+1]   =   p1*Amax*(float_type)exp(-T*T)*(exp(j*carrier_phase)+noiselevel*(frand0()+j*frand0()))+A0*noiselevel_uni*(frand0()+j*frand0()); 
     }

     fftw(FFT_FWPLAN_T, 2, (fftwt_complex*)FIELD_REAL_SMALL, 2, 1, (fftwt_complex*)(NL_SMALLBUFFER1), 2,1); 
	 for (int cw=0; cw<2*N_T; cw++)
	 {
	  float_type phase = (R*R*real(WAVENUMBER[cw])/2.0/f);
	  //float_type phase = (R*R*real(WAVENUMBER0[cw%2])/2.0/f);
	  FIELD[ofs0+cw]   = (NL_SMALLBUFFER1[cw]*exp(j*phase));
	  //FIELD[ofs0+cw]   = NL_SMALLBUFFER1[cw];
	 } 
        }
        if (ISMASTER) printf("Done.");
}
