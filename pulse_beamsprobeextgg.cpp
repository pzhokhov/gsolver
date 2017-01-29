#include "solver.h"
#include "pulseshapes.h"

void create_beamsprobeextgg(FILE* fid)
{
    f_complex j = f_complex(0,1);
    char pulseshape[200] = ""; 
    load_namedstringn(fid, "PULSE_SHAPE", pulseshape, 200);
    float_type noiselevel = load_namedfloat(fid, "NOISE_LEVEL", true, 0);
    float_type noiselevel_uni = load_namedfloat(fid, "NOISE_LEVEL_UNIFORM", true, 0);
	float_type spatialnoiselevel = load_namedfloat(fid, "NOISE_LEVEL_SPATIAL", true, 0);
    
   float_type probepower_ratio = load_namedfloat(fid, "PROBE_POWERRATIO",true,1);
    float_type probefield_ratio = sqrt(probepower_ratio);

    float_type Pcr   = 2.0*M_PI*LIGHT_VELOCITY/real(WAVENUMBER0)/OMEGA0/NONLIN_REFRINDEX*INTENSITY_DENOM;
    
          
    int Nbeams = atoi(pulseshape); 
    float_type* theta     = (float_type*)malloc_ch(Nbeams*sizeof(float_type)); 
    float_type* phi       = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
    float_type* d_fwhm    = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
    float_type* tau_fwhm  = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
    float_type* E         = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
    float_type* foc       = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
    float_type* crshift   = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
    float_type* omega     = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
    float_type* delay     = (float_type*)malloc_ch(Nbeams*sizeof(float_type));	
    float_type* cephase   = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
	float_type* C         = (float_type*)malloc_ch(Nbeams*sizeof(float_type)); 
	float_type* Cw		  = (float_type*)malloc_ch(Nbeams*sizeof(float_type));
    float_type* maskA     = (float_type*)malloc_ch(Nbeams*sizeof(float_type)); 
    float_type* maskNx    = (float_type*)malloc_ch(Nbeams*sizeof(float_type)); 
    float_type* maskNy    = (float_type*)malloc_ch(Nbeams*sizeof(float_type)); 



    float_type fmean = 0;   // mean focusing distance. It is used for inverse diffraction part. 
    for (int i=0; i<Nbeams; i++)
    {
     theta[i]    = load_namednumfloat(fid,"THETA",i);
     phi[i]      = load_namednumfloat(fid,"PHI",i,true,2*M_PI/Nbeams*i);
     d_fwhm[i]   = load_namednumfloat(fid,"FOCUS_WAIST_LINEAR_FWHM", i, true, -1);
     tau_fwhm[i] = load_namednumfloat(fid, "DURATION_FWHM", i);
     E[i]        = load_namednumfloat(fid, "ENERGY",i); 
     crshift[i]  = load_namednumfloat(fid, "CROSSING_SHIFT",i, true, 0); 
     foc[i]      = load_namednumfloat(fid, "FOCUSING_DISTANCE",i);
     omega[i]    = 2*M_PI*LIGHT_VELOCITY/load_namednumfloat(fid, "LAMBDA_V",i); 
     delay[i]    = load_namednumfloat(fid, "DELAY",i,true, 0);
     cephase[i]  = load_namednumfloat(fid, "CE_PHASE",i,true, 2.0*M_PI*frand_common());
	 Cw[i]       = load_namednumfloat(fid, "CHIRP_SPECTRAL", i, true, 0);
	 C[i]        = load_namednumfloat(fid, "CHIRP_SI", i, true, 0);
	 C[i]        = load_namednumfloat(fid, "CHIRP",    i, true, C[i]*tau_fwhm[i]*tau_fwhm[i]);
	     
	maskA[i] = load_namednumfloat(fid, "MASK_A", i, true, 0); 
	maskNx[i] = load_namednumfloat(fid, "MASK_NX", i, true, 1); 
	maskNy[i] = load_namednumfloat(fid, "MASK_NY", i, true, 1);        


     if (d_fwhm[i] == -1)
     {
      // beam waist at focus point is not defined, let us load initial beam diameter
      float_type a = load_namednumfloat(fid, "DIAMETER_FWHM", i)/sqrt(log(4.0));
      float_type zR = real(WAVENUMBER0)*a*a/2.0;
      float_type f=foc[i]; if (zR > 2*f) f = zR*zR  /2/foc[i]*(1-sqrt(1-4*foc[i]*foc[i]/zR/zR));
      d_fwhm[i] = a*sqrt(log(4.0))*f/zR/sqrt(1+(f*f/zR/zR)); // calculate linear beam waist at focus point from initial diameter and focusing distance
     }
     fmean += foc[i];
    }
    fmean /= Nbeams;
    for (int i=0; i<Nbeams; i++) foc[i]-=fmean;

    if (ISMASTER)
    {
     printf("\nCreating %d gaussian beams, last(probe) being %f times weaker...\n distance from starting plane to focus point %e", Nbeams, probepower_ratio, fmean);
     for (int i=0; i < Nbeams; i++) printf("\n %d)theta=%f, phi=%f, d_fwhm=%e, tau_fwhm=%e, E=%e, foc=%e, crshift=%e, cephase=%g, C=%g, Cw=%g, delay=%g",i,theta[i], phi[i], d_fwhm[i], tau_fwhm[i], E[i], foc[i], crshift[i], cephase[i], C[i], Cw[i], delay[i]);
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
      float_type s = sin(theta[i]), c = cos(theta[i]);
      float_type A0 = sqrt(E[i]/d_fwhm[i]/d_fwhm[i]/tau_fwhm[i]*sqrt(log(16.0)/M_PI)*log(16.0)/M_PI);              
      float_type tau = tau_fwhm[i]/sqrt(log(4.0));
      float_type ro  = d_fwhm[i]/sqrt(log(4.0));
      float_type zR    = real(WAVENUMBER0)*ro*ro/2;
      float_type Pin   = E[i]/tau_fwhm[i]*sqrt(log(16.0)/M_PI); 
      float_type delta = 0;
      float_type f = foc[i]+fmean;
      if (Pin > Pcr) delta = 0; //f*f/(f+0.367*zR*(1+f*f/zR/zR)/sqrt((sqrt(Pin/Pcr)-0.852)*(sqrt(Pin/Pcr)-0.852)-0.0219));
      float_type D     = -s*delta + crshift[i]*s/c;
      float_type alldelay = delay[i];//(f-delta+crshift[i]*cos(theta[i]))/GROUP_VELOCITY*(cos(theta[i])-cos(theta[0]));

      float_type xlocal = x - D*cos(phi[i]), ylocal = y - D*sin(phi[i]);       
      float_type X = ((xlocal*cos(phi[i])+ylocal*sin(phi[i])));//*cos(theta[i]));
      float_type Y = (-xlocal*sin(phi[i])+ylocal*cos(phi[i]));
      float_type R = sqrt((X*X+Y*Y))/ro;
//      float_type Amax = A0*exp(-R*R);
			
//      foc[i] = 0;
      float_type X_ =  X*c - foc[i]*s;
      float_type Z_ = -X*s - foc[i]*c;
	  f_complex A0_ = A0;//*frand1c(spatialnoiselevel); 
      for (int ct=0; ct<N_T; ct++) 
      {
       float_type t = TMIN + TSTEP*ct;
       float_type T = ((-Z_-foc[i]/cos(theta[0]))/GROUP_VELOCITY +  t - alldelay)/tau;
       
       float_type carrier_phase = -(omega[i]-OMEGA0)*t+cephase[i] - ((OMEGA0 == OMEGA[0])?(0):(OMEGA0*t)) - C[i]*T*T - OMEGA0*delay[i];
					
       FIELD_REAL_SMALL[ct] = A0_*(float_type)exp(-T*T)*(exp(j*carrier_phase)/*+noiselevel*(frand0()+j*frand0())*/); 
       if (i==Nbeams-1) FIELD_REAL_SMALL[ct]  *=probefield_ratio;
       // FIELD_REAL_SMALL[ct]   += A0_*noiselevel_uni*(frand0()+j*frand0());
      }
      fftwt_execute_dft(FFT_FWPLAN_T, (fftwt_complex*)FIELD_REAL_SMALL, (fftwt_complex*)NL_SMALLBUFFER1); 
      for (int cw=0; cw<N_T; cw++)
      {
//	float_type phase = (R*R*ro*ro*real(WAVENUMBER[cw])/2/f)+X*tan(theta[i])*real(WAVENUMBER[cw]);
//      float_type phase = ((x*x + y*y)*real(WAVENUMBER[cw])/2/f);
        //*(1+maskA[i]*sin(2*M_PI*X/ro/maskNx[i])*sin(2*M_PI*X/ro/maskNx[i]))

        f_complex zR_ = WAVENUMBER[cw]*ro*ro/(float_type)2.0;
        f_complex cphase = -(X_*X_ + Y*Y)/(ro*ro*((float_type)1.0-j*Z_/zR)) - j*WAVENUMBER[cw]*Z_ - j*OMEGA[cw]*foc[i]/GROUP_VELOCITY/(float_type)cos(theta[0]) - j*Cw[i]*(OMEGA[cw]-OMEGA0)*(OMEGA[cw]-OMEGA0)/((float_type)2.0);
	f_complex mult   = (float_type)1.0/((float_type)1.0-j*Z_/zR_);
	FIELD[ofs0+cw] += (NL_SMALLBUFFER1[cw]*exp(cphase)*mult);
      }
     }
    }
    
    //Inverse diffraction part. 
    GROUP_VELOCITY *= cos(theta[0]);
     
    MPI_Barrier(MPI_COMM_WORLD);
    fftwt_mpi_execute_dft(FFT_FWPLAN_XY, (fftwt_complex*)FIELD,      (fftwt_complex*)BIGBUFFER2); 
                    
    float_type kxstep = 2*M_PI/(XMAX-XMIN), kystep = 2*M_PI/(YMAX-YMIN);
    float_type maxI = 0;

    for (int nx = 0; nx<MY_NX_FT; nx++)
    {
     int nx_g = nx+MY_NXstart_FT;   
     float_type kx=0; if (nx_g<=N_X/2) kx=kxstep*nx_g; else kx=-kxstep*(N_X-nx_g);

     for (int ny = 0; ny<N_Y;     ny++)
     {           
      float_type ky=0; if (ny<=N_Y/2) ky=kystep*ny; else ky=-kystep*(N_Y-ny); 
      float_type kt2 = kx*kx+ky*ky;
  	
      for (int nw = 0; nw<N_T;     nw++)
      {
       int ofs = (nw + N_T*(ny+N_Y*nx)); 
       if (OMEGA[nw] < 0 || kt2 > (MAX_KT2*real(WAVENUMBER[nw]*WAVENUMBER[nw])))     {FIELD[ofs]=0.0;continue;}
             
       //kz is like in linear propagator part, witout imaginary part however
	   f_complex k1 = WAVENUMBER0 + OMEGA[nw]/GROUP_VELOCITY;
       float_type kz = real(sqrt(k1*k1 - kt2)); 
	     
       FIELD[ofs] *= exp(-j*fmean*(OMEGA[nw]/GROUP_VELOCITY - kz));
	   maxI = max(maxI, abs2(FIELD[ofs]));
      }
     }
    }   
    MPI_Barrier(MPI_COMM_WORLD);
    fftwt_execute(FFT_BWPLAN_XY); 
    for (size_t i=0; i<N_T*N_X*MY_NY; i++) FIELD[i] /= (N_X*N_Y);  //fftw N-normalization*/

	//Noise addition section
	fftwt_execute(FFT_ALLBWPLAN_T); fftwt_Nnormalize(MY_NX*MY_NY, BIGBUFFER1);
    maxI = 0;
	for (int ny=0; ny<MY_NY; ny++) for (int nx=0; nx<N_X; nx++) 
	{
	 long int ofs0 = N_T*(nx+N_X*ny);
	 for (int nt=0; nt<N_T; nt++)
	 {
		long int ofs = nt+ofs0;
		maxI = max(maxI, abs2(BIGBUFFER1[ofs]));
		BIGBUFFER1[ofs] *= (frand1c(noiselevel));
	 }
	 f_complex r = frand1c(spatialnoiselevel);
	 for (int nt=0; nt<N_T; nt++) BIGBUFFER1[nt+ofs0] *= r;
	}

	fftwt_plan fft_allfwplan_t = fftwt_plan_many_dft(1, &N_T, MY_NX*MY_NY, (fftwt_complex*)BIGBUFFER1, NULL, 1, N_T, (fftwt_complex*)FIELD, NULL, 1, N_T, FFTW_FORWARD, FFTW_FLAG); 
	fftwt_execute(fft_allfwplan_t);

    free(theta); free(phi); free(d_fwhm); free(tau_fwhm); free(E);
    free(foc); free(crshift);
    free(omega); free(delay);
	free(C); free(Cw);
    free(cephase);
    free(maskA); free(maskNx); free(maskNy); 
}	

