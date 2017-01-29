#include "solver.h"
#include "pulseshapes.h"

void load_spectrum(f_complex* output);

void create_customspectrumg(FILE* fid)
{
        f_complex j = f_complex(0,1);
    	float_type d_fwhm   = load_namedfloat(fid, "DIAMETER_FWHM", true, sqrt(log(16.0))*load_namedfloat(fid, "RADIUS_GAUSS", true, 0));
		if (d_fwhm == 0) throw "Neither of parameters specifying beam radius is given!";

    	float_type E        = load_namedfloat(fid, "ENERGY");
    	float_type noiselevel = load_namedfloat(fid, "NOISE_LEVEL", true, 0);
		float_type noiselevel_uni = load_namedfloat(fid, "NOISE_LEVEL_UNIFORM",true,0);
		float_type spatialnoiselevel = load_namedfloat(fid, "NOISE_LEVEL_SPATIAL", true, 0);

		load_spectrum(NL_SMALLBUFFER1);
	

	float_type phase_flip        = load_namedfloat(fid, "PHASE_FLIP", true, 0)*M_PI; 
	float_type phase_dip	     = load_namedfloat(fid, "PHASE_DIP",  true, 0);
  	float_type phase_flip_l      = load_namedfloat(fid, "PHASE_FLIP_LAMBDA", true, 1); 
	float_type phase_flip_width  = load_namedfloat(fid, "PHASE_FLIP_WIDTH",  true, 0); 

	float_type phase_flip_maxw = 2*M_PI*LIGHT_VELOCITY/(phase_flip_l - phase_flip_width/2.0); 
	float_type phase_flip_minw = 2*M_PI*LIGHT_VELOCITY/(phase_flip_l + phase_flip_width/2.0); 


	
        if (ISMASTER) printf("\nCreating gaussian-beam custom spectrum pulse...");

		float_type ro  = d_fwhm/sqrt(log(4.0));
        float_type A0 = sqrt(E/d_fwhm/d_fwhm*log(16.0)/M_PI*N_T/(TMAX-TMIN));              
 
        
        float_type zR  = real(WAVENUMBER0)*ro*ro/2;
            
		float_type f  = load_namedfloat(fid, "FOCUSING_DISTANCE", true, 0);
        if (zR > 2*f) f = zR*zR  /2/f*(1-sqrt(1-4*f*f/zR/zR)); //inital wavefront curvature radius calculated from distance from linear focus
		f        = load_namedfloat(fid, "WAVEFRONT_RADIUS", true, f);	            

        if (ISMASTER)
        {
         printf("\n   A0  = %e", A0);
         printf("\n   ro  = %e", ro);
		 printf("\n   f   = %e", f);
        }
                        
	for (int cy=0; cy<MY_NY; cy++)
	for (int cx=0; cx<MY_NX; cx++)
	{
	 int ofs0 = N_T*(cx+MY_NX*cy);                
     float_type R; 

#ifdef _UNIAXIAL
	 R = HT_PLAN->x_n(cx+MY_NXstart)*XMAX;
#else
	 float_type x = XMIN+cx*XSTEP, y = YMIN+(cy+MY_NYstart)*YSTEP;
	 R = sqrt(x*x+y*y);
#endif
	 f_complex A0_ = A0*frand1c(spatialnoiselevel); 
     for (int cw=0; cw<N_T; cw++) 
     {	
	  	float_type phase = 0; 
          	if (phase_flip != 0 && phase_flip_minw < OMEGA[cw] && OMEGA[cw] < phase_flip_maxw) phase += phase_flip; 
		FIELD[ofs0+cw] = exp(-R*R/ro/ro + j*R*R*real(WAVENUMBER[cw])/(float_type)2.0/f -j*phase)*A0_*(NL_SMALLBUFFER1[cw] + noiselevel*frand0c()) + A0*noiselevel_uni*frand0c();
		if (phase_dip != 0  && phase_flip_minw < OMEGA[cw] && OMEGA[cw] < phase_flip_maxw) FIELD[ofs0+cw] *= (1.0-phase_dip); 
     }
	}
	if (ISMASTER) {printf("Done."); fflush(stdout);}
}



void load_spectrum(f_complex* output)
{
	char spectrumfilename[512]; 
	char spectrumfiletype[128];

	f_complex j = f_complex(0,1);

	load_namedstringn(SC_FID, "INPUT_SPECTRUM_FILE", spectrumfilename, 512);
	load_namedstringn(SC_FID, "INPUT_SPECTRUM_TYPE", spectrumfiletype,  128); 

	float_type Cw = load_namedfloat(SC_FID, "CHIRP_SPECTRAL", true , 0);

	load_spectrum_fromfile(output,OMEGA,N_T, spectrumfilename, spectrumfiletype);

	for (int nw=0; nw<N_T; nw++) 
 	{
	  output[nw] *= exp(-j*((OMEGA[nw]-OMEGA0)*(OMEGA[nw]-OMEGA0)*Cw*(float_type)0.5+OMEGA[nw]*TMIN));	 
 	}
}


