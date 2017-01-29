#include "solver.h"
#include "pulseshapes.h"

void create_beamsprobegg(FILE* fid)
{
    f_complex j = f_complex(0,1);
    char pulseshape[200] = ""; 
    load_namedstringn(fid, "PULSE_SHAPE", pulseshape, 200);
    float_type tau_fwhm = load_namedfloat(fid, "DURATION_FWHM"); 
    float_type d_fwhm   = load_namedfloat(fid, "DIAMETER_FWHM");
    float_type E        = load_namedfloat(fid, "ENERGY");
    float_type f        = load_namedfloat(fid, "FOCUSING_DISTANCE");
    float_type noiselevel = load_namedfloat(fid, "NOISE_LEVEL", true, 0);
    float_type noiselevel_uni = load_namedfloat(fid, "NOISE_LEVEL_UNIFORM", true, 0);
    int Nbeams = atoi(pulseshape); 
    float_type theta = atof(pulseshape+strlen(SHAPE_BEAMSPROBEGG)+1);
    
    float_type A0 = sqrt(E/d_fwhm/d_fwhm/tau_fwhm*sqrt(log(16.0)/M_PI)*log(16.0)/M_PI/cos(theta/2.0));              
    float_type tau = tau_fwhm/sqrt(log(4.0));
    float_type ro  = d_fwhm/sqrt(log(4.0));
    float_type probepower_ratio = load_namedfloat(fid, "PROBE_POWERRATIO");
    float_type probefield_ratio = sqrt(probepower_ratio);
            
            
    float_type zR    = real(WAVENUMBER0)*ro*ro/2;
    float_type Pcr   = 2.0*M_PI*LIGHT_VELOCITY/real(WAVENUMBER0)/OMEGA0/NONLIN_REFRINDEX;
    float_type Pin   = E/tau_fwhm*sqrt(log(16.0)/M_PI);
    float_type delta = 0;
    if (Pin > Pcr) delta = f*f/(f+0.367*zR/sqrt((sqrt(Pin/Pcr)-0.852)*(sqrt(Pin/Pcr)-0.852)-0.0219));

    float_type D     = 2.0*sin(theta/2.0)*(f-delta) + load_namedfloat(fid, "CROSSING_SHIFT", true, 0)*2.0*tan(theta/2.0);
   
    if (zR > 2*f) f = zR*zR  /2/f*(1-sqrt(1-4*f*f/zR/zR));
   	
    if (ISMASTER)
    {
     printf("\nCreating %d gaussian beams, last(probe) being %f times weaker...", Nbeams, probepower_ratio);
     printf("\n A0 = %e", A0);
     printf("\n tau = %e", tau);
     printf("\n ro = %e",  ro);
     printf("\n theta = %e", theta);
     printf("\n D  = %e", D);
     printf("\n Pcr = %e", Pcr);
     fflush(stdout);
    }   
            
            	
    for (int cy=0; cy<MY_NY; cy++)
    for (int cx=0; cx<N_X; cx++)
    {
     float_type x = XMIN+cx*XSTEP, y = YMIN+(cy+MY_NYstart)*YSTEP;
     int ofs0 = N_T*(cx+N_X*cy);
     fill_zeros(FIELD+ofs0, N_T); 
     for (int i=0; i<Nbeams; i++)
     {
      float_type phi = 2*M_PI*i/Nbeams;
      float_type xlocal = x - D/2*cos(phi), ylocal = y - D/2*sin(phi);       
      float_type X = ((xlocal*cos(phi)+ylocal*sin(phi))*cos(theta/2));
      float_type Y = (-xlocal*sin(phi)+ylocal*cos(phi));
      float_type R = sqrt((X*X+Y*Y))/ro;
      float_type Amax = A0*exp(-R*R);
				
      for (int ct=0; ct<N_T; ct++) 
      {
       float_type t = TMIN + TSTEP*ct;
       float_type T = ((xlocal*cos(phi)+ylocal*sin(phi))/GROUP_VELOCITY*sin(theta/2)   +  t)/tau;
			
       float_type carrier_phase = (OMEGA0 == OMEGA[0])?(0):(-OMEGA0*(TMIN+ct*(TMIN-TMAX)/N_T));
					
       FIELD_REAL_SMALL[ct] = Amax*(float_type)exp(-T*T)*(exp(j*carrier_phase)+noiselevel*(frand0()+j*frand0())); 
       if (i==Nbeams-1) FIELD_REAL_SMALL[ct] *= probefield_ratio;
       FIELD_REAL_SMALL[ct] += A0*noiselevel_uni*(frand0()+j*frand0());
      }

      fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)FIELD_REAL_SMALL, (fftwt_complex*)NL_SMALLBUFFER1); 
      for (int cw=0; cw<N_T; cw++)
      {
       float_type phase = (R*R*ro*ro*real(WAVENUMBER[cw])/2/f)+X*tan(theta/2)*real(WAVENUMBER[cw]);
       FIELD[ofs0+cw] += (NL_SMALLBUFFER1[cw]*exp(j*phase));
      }
     }
    }
    GROUP_VELOCITY *= cos(theta/2);
}	

