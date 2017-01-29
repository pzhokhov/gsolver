#include "dumps.h"
#include "solver.h"

dump_type_empty_class :: ~dump_type_empty_class()
{
	fclose(dumpfid);
}


void dump_type_empty_class :: init_dumpfile(FILE* optfid, int i)
{
	char filenamebuf[512];
	char filenamebuf1[512];
	char modnamebuf[100];
	load_namednumstringn(optfid, "DUMP_FILENAME", filenamebuf1, i, 200, false);
	load_namedstringn(optfid, "DUMP_FILENAME_MOD", modnamebuf, 100, true, "");
	sprintf(filenamebuf, filenamebuf1, modnamebuf);
	if (ISMASTER) 
	 printf("\nInitializing dump file %s (tN=%d, xN=%d, yN=%d, zN=%d)...", filenamebuf, tN, xN, yN, zN);

	init_dumpfile(filenamebuf);
	if (ISMASTER) {printf("Done."); }
}

void dump_type_empty_class :: init_dumpfile(const char* filename)
{
	const char dumpstandard_name[] = "UDF1.0";
	int float_size = sizeof(float_type);

	if (ISMASTER)
	{
	 dumpfid = fopen(filename, "wb");
	 if (!dumpfid) 
	 {
		printf("\nvoid dump_type_empty_class :: init_dumpfile: Unable to open file %s for writing.",filename); 
		perror("Unable to open file!"); fflush(stdout);
		throw "void dump_type_empty_class :: init_dumpfile : Unable to open file!";
	 }
	 fwrite(dumpstandard_name, sizeof(char), strlen(dumpstandard_name)+1, dumpfid);
	 fwrite(&float_size, sizeof(int), 1, dumpfid);
	 fwrite(&iscomplex,  sizeof(int), 1, dumpfid);
	 fwrite(&tN,         sizeof(int), 1, dumpfid);
	 fwrite(&xN,         sizeof(int), 1, dumpfid);
	 fwrite(&yN,         sizeof(int), 1, dumpfid);
	 fwrite(&zN,         sizeof(int), 1, dumpfid);
	 fclose(dumpfid);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	dumpfid=fopen(filename, "r+b"); 

	header_ofs = sizeof(char)*(strlen(dumpstandard_name)+1) + 6*sizeof(int);
}

void dump_type_empty_class :: dump()
{
	dump_noflush();
	fflush(dumpfid);
}

void dump_type_empty_class :: dump_write(void* buf)
{
    long long int ofs = header_ofs; 
	if (iscomplex) ofs += (long long int)fileofs()*sizeof(f_complex); else ofs += (long long int)fileofs()*sizeof(float_type);
#ifdef _WIN32
	_fseeki64(dumpfid, ofs, SEEK_SET);
#else
	fseeko(dumpfid, ofs, SEEK_SET);
#endif

	if (iscomplex) fwrite(buf, sizeof(f_complex) , piece_size, dumpfid);
	else           fwrite(buf, sizeof(float_type), piece_size, dumpfid);
}




dump_type_full_class :: dump_type_full_class(FILE* fid, int num)
{
	nz_denom = load_namednumint(fid, "NZ_DENOM", num, true, 1);
	xN = N_X; yN=N_Y; zN=N_Z/nz_denom; tN=N_T;  
	piece_size = MY_NX*tN*MY_NY;
	iscomplex = 1; 
	init_dumpfile(fid, num);
}


long int dump_type_full_class :: fileofs()
{
	return tN*(xN*yN*(n_Z/nz_denom) + MY_NXstart + N_X*MY_NYstart);
}

void dump_type_full_class :: dump_noflush()
{
	if (n_Z%nz_denom) return;
	fftwt_execute(FFT_ALLBWPLAN_T); fftwt_Nnormalize(MY_NX*MY_NY, BIGBUFFER1);
	dump_write(BIGBUFFER1); 
}


dump_type_maxI_class :: dump_type_maxI_class(FILE* fid, int num)
{
	xN = N_X; yN=N_Y; zN=N_Z; tN=1; 
	piece_size = MY_NX*MY_NY;
	iscomplex = 0; 
	init_dumpfile(fid, num);

	wfilter = (f_complex*)malloc_ch(sizeof(f_complex)*N_T); 
	char filtertype[100];  load_namednumstringn(fid, "FILTER_TYPE", filtertype, num, 100, "true", "NO");
	create_wfilter(wfilter, filtertype);
}

dump_type_maxI_class :: ~dump_type_maxI_class() 
{
	free(wfilter);
}

long int dump_type_maxI_class :: fileofs()
{
	return xN*(yN*n_Z+MY_NYstart) + MY_NXstart;
}

void dump_type_maxI_class :: dump_noflush()
{
	float_type* buf = (float_type*)BIGBUFFER1;
	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<MY_NX; nx++)
	{
		int ofs = N_T*(nx+MY_NX*ny);
		for(int nw=0; nw<N_T; nw++) NL_SMALLBUFFER1[nw]=FIELD[ofs+nw]*wfilter[nw]/(float_type)N_T;
		fftwt_execute_dft(FFT_BWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)(BIGBUFFER1+ofs)); 
	}

	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<MY_NX; nx++)
	{
		float_type maxI=0;
		for(int nt=0; nt<N_T; nt++) 
		{
			int ofs = nt + N_T*(nx+MY_NX*ny); maxI = max(maxI, abs2(BIGBUFFER1[ofs]));
		}
		buf[nx+MY_NX*ny] = maxI;
	}
	dump_write(buf);
}

void dump_type_flux_class :: dump_noflush()
{
	// TODO: Compare results for frequency- and time-domain integration
	float_type* buf = (float_type*)BIGBUFFER1;
	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<MY_NX; nx++)
	{
		int ofs = N_T*(nx+MY_NX*ny);
		for(int nw=0; nw<N_T; nw++) NL_SMALLBUFFER1[nw]=FIELD[ofs+nw]*wfilter[nw]/(float_type)N_T;
		fftwt_execute_dft(FFT_BWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)(BIGBUFFER1+ofs)); 
	}

	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<MY_NX; nx++)
	{
		float_type f=0;
		for(int nt=0; nt<N_T; nt++) 
		{
			int ofs = nt + N_T*(nx+MY_NX*ny); f+=abs2(BIGBUFFER1[ofs]);
		}
		buf[nx+MY_NX*ny] = f*TSTEP;
	}
	dump_write(buf);
}

void dump_type_duration_class :: dump_noflush()
{
	float_type* buf = (float_type*)BIGBUFFER1;
	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<MY_NX; nx++)
	{
		int ofs = N_T*(nx+MY_NX*ny);
		for(int nw=0; nw<N_T; nw++) NL_SMALLBUFFER1[nw]=FIELD[ofs+nw]*wfilter[nw];
		fftwt_execute_dft(FFT_BWPLAN_T, (fftwt_complex*)NL_SMALLBUFFER1, (fftwt_complex*)(BIGBUFFER1+ofs)); 
	}
    fftwt_Nnormalize(MY_NX*MY_NY, BIGBUFFER1);

	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<MY_NX; nx++)
	{
		float_type f=0, ft=0, ft2=0;
		for(int nt=0; nt<N_T; nt++) 
		{
			float_type t = TMIN+TSTEP*nt;
			int ofs = nt + N_T*(nx+MY_NX*ny);
			float_type I =abs2(BIGBUFFER1[ofs]);
			f   += I;
			ft  += I*t;
			ft2 += I*t*t; 
		}
		ft  /= f;
		ft2 /= f; 

		buf[nx+MY_NX*ny]=sqrt(ft2-ft*ft);
	}
	dump_write(buf);
}

dump_type_fluxk_class :: dump_type_fluxk_class(FILE* fid, int num)
{
	nz_denom = load_namednumint(fid, "NZ_DENOM", num, true, 1);
	xN = N_Y; yN=N_X; zN=N_Z/nz_denom; tN=1; 
	piece_size = tN*N_Y*MY_NX_FT;
	iscomplex = 0; 
	init_dumpfile(fid, num);

	wfilter = (f_complex*)malloc_ch(sizeof(f_complex)*N_T); 
	char filtertype[100];  load_namednumstringn(fid, "FILTER_TYPE", filtertype, num, 100, "true", "NO");
	create_wfilter(wfilter, filtertype);
}

long int  dump_type_fluxk_class :: fileofs()
{
	return tN*xN*(yN*n_Z/nz_denom+MY_NXstart_FT);
}


void dump_type_fluxk_class :: dump_noflush()
{
	if (n_Z%nz_denom) return;
	float_type* buf = (float_type*)BIGBUFFER1;
	for(int ny=0; ny<N_Y; ny++) for(int nx=0; nx<MY_NX_FT; nx++)
	{
		float_type F = 0; 
		int ofs = N_T*(ny+N_Y*nx);
		for(int nw=0; nw<N_T; nw++) F += abs2(BIGBUFFER2[ofs+nw]*wfilter[nw]);
		buf[ny+N_Y*nx] = F;
	}
	dump_write(buf);
}

void dump_type_ysection_class :: init_ny(FILE* optfid, int num)
{
	ny = load_namednumint(optfid, "DUMP_NY", num, true, N_Y/2);
	ny -= MY_NYstart; if (0 <= ny && ny < MY_NY) mydump=true; else mydump=false;
}

dump_type_ysection_class :: dump_type_ysection_class(FILE* optfid, int num)
{
	tN=N_T; xN = N_X; yN=1; zN=N_Z;
	piece_size=tN*MY_NX; iscomplex=1;
	init_dumpfile(optfid, num);
	init_ny(optfid, num);
	if (mydump) fft_bwplan_sectiont = fftwt_plan_many_dft(1, &N_T, N_X, (fftwt_complex*)(FIELD+N_T*MY_NX*ny), NULL, 1, N_T,  (fftwt_complex*)BIGBUFFER1, NULL, 1, N_T, FFTW_BACKWARD, FFTW_FLAG); 
}

long int dump_type_ysection_class :: fileofs()
{
	return tN*xN*n_Z + MY_NXstart;
}

void dump_type_ysection_class :: dump_noflush()
{
	MPI_Barrier(MPI_COMM_WORLD); 
	if (mydump) 
	{
		fftwt_execute(fft_bwplan_sectiont); fftwt_Nnormalize(MY_NX, BIGBUFFER1);
		dump_write(BIGBUFFER1);
	}
	MPI_Barrier(MPI_COMM_WORLD); 
}

dump_type_ysectionk_class :: dump_type_ysectionk_class(FILE* optfid, int num) : dump_type_ysection_class() 
{
	tN=N_T; xN = N_X; yN=1; zN=N_Z;
	piece_size=tN*MY_NX_FT; iscomplex=1;
	init_dumpfile(optfid, num);
	init_ny(optfid, num);
}


void dump_type_ysectionk_class :: init_ny(FILE* optfid, int num)
{
	ny = load_namednumint(optfid, "DUMP_NY", num, true, 0);
}

long int dump_type_ysectionk_class :: fileofs()
{
	return (long int)tN*(long int)(xN*n_Z+MY_NXstart_FT);
}

void dump_type_ysectionk_class :: dump_noflush()
{
	for (size_t nx=0; nx<MY_NX_FT; nx++) 
		for (size_t nw=0; nw<N_T;      nw++) 
			BIGBUFFER1[nw+N_T*nx] = BIGBUFFER2[nw+N_T*(ny+N_Y*nx)];
	dump_write(BIGBUFFER1);
}



dump_type_ysection_maxI_class :: dump_type_ysection_maxI_class(FILE* optfid, int num)
{
	tN=1; xN = N_X; yN=1; zN=N_Z;
	piece_size=xN; iscomplex=0;
	init_dumpfile(optfid, num);
	init_ny(optfid, num);
	fft_bwplan_sectiont = fftwt_plan_many_dft(1, &N_T, N_X, (fftwt_complex*)FIELD+N_T*MY_NX*ny, NULL, 1, N_T,  (fftwt_complex*)BIGBUFFER1, NULL, 1, N_T, FFTW_BACKWARD, FFTW_FLAG); 

	wfilter = (f_complex*)malloc_ch(sizeof(f_complex)*N_T); 
	char filtertype[100];  load_namednumstringn(optfid, "FILTER_TYPE", filtertype, num, 100, "true", "NO");
	create_wfilter(wfilter, filtertype);
}

dump_type_ysection_maxI_class :: ~dump_type_ysection_maxI_class()
{
	free(wfilter);
}

void dump_type_ysection_maxI_class :: dump_noflush()
{
	if (mydump) 
	{
		float_type* buf = (float_type*)BIGBUFFER1;
		fftwt_execute(fft_bwplan_sectiont); fftwt_Nnormalize(N_X, BIGBUFFER1);
		for (int nx=0; nx<N_X; nx++) 
		{
			float_type maxI = 0;
			for (int nt=0; nt<N_T; nt++) maxI=max(maxI, abs2(BIGBUFFER1[nt+N_T*nx]));
			buf[nx] = maxI;
		}
		dump_write(buf);
	}
}



dump_type_youngy_class :: dump_type_youngy_class(FILE* optfid, int num)
{
	ny = N_Y/2;
	tN=1; xN = N_Y/2; yN=N_X; zN=N_Z;
	piece_size=MY_NX_FT*N_Y/2; iscomplex=1;
	init_dumpfile(optfid, num);
	

	wfilter = (f_complex*)malloc_ch(sizeof(f_complex)*N_T); 
	char filtertype[100];  load_namednumstringn(optfid, "FILTER_TYPE", filtertype, num, 100, "true", "NO");
	create_wfilter(wfilter, filtertype);
}

dump_type_youngy_class :: ~dump_type_youngy_class()
{
	free(wfilter);
}

void dump_type_youngy_class :: dump_noflush()
{
	 for (size_t nx=0;  nx<MY_NX_FT; nx++)
	 for (size_t ny_=0; ny_<N_Y/2; ny_++)
	 {
      f_complex mA12 = 0;
      float_type mA1 = 0, mA2 = 0;
	  for (size_t nw=0;  nw<N_T; nw++)
	  {
		 size_t ofs1 = nw+N_T*(      ny_+N_Y*nx);
		 size_t ofs2 = nw+N_T*(N_Y-1-ny_+N_Y*nx); 
		 f_complex A1 = BIGBUFFER2[ofs1]*wfilter[nw], A2=BIGBUFFER2[ofs2]*wfilter[nw];
		 mA1 += abs2(A1);
		 mA2 += abs2(A2);
		 mA12 += A1*conj(A2);
	  }
	  BIGBUFFER1[ny_+N_Y/2*nx]= ((float_type)2.0)*(mA12)/(mA1+mA2);
      //BIGBUFFER1[ny_+N_Y/2*nx]=mA1; 
	 }
	 dump_write(BIGBUFFER1);
}


long int dump_type_youngy_class :: fileofs()
{
  return tN*xN*(yN*n_Z + MY_NXstart_FT);
}


void dump_type_ysection_flux_class :: dump_noflush()
{
	if (mydump) 
	{
		float_type* buf = (float_type*)BIGBUFFER1;
		fftwt_execute(fft_bwplan_sectiont); fftwt_Nnormalize(N_X, BIGBUFFER1);
		for (int nx=0; nx<N_X; nx++) 
		{
			float_type f = 0;
			for (int nt=0; nt<N_T; nt++) f+=abs2(BIGBUFFER1[nt+N_T*nx]);
			buf[nx] = f*TSTEP;
		}
		dump_write(buf);
	}
}


dump_type_field_axis_class :: dump_type_field_axis_class(FILE* optfid, int num)
{
	tN = N_T; xN = 1; yN = 1; zN = N_Z;
	piece_size=N_T; iscomplex=1;
	init_dumpfile(optfid, num);
}


long int dump_type_field_axis_class :: fileofs()
{
	return tN*n_Z;
}

void dump_type_field_axis_class :: dump_noflush()
{
#ifdef _UNIAXIAL
	if (!ISMASTER) return;
	int nx = 0; 
	int ny = 0; 
#else
	int ny = N_Y/2 - MY_NYstart; 
	if ((ny<0) || (MY_NY<=ny)) return;
	int nx = N_X/2;
#endif
	fftwt_execute_dft(FFT_BWPLAN_T, (fftwt_complex*)FIELD+N_T*(nx+N_X*ny), (fftwt_complex*)BIGBUFFER1); fftwt_Nnormalize(1, BIGBUFFER1);
	dump_write(BIGBUFFER1);
}


dump_type_average_spectrum_class :: dump_type_average_spectrum_class(FILE* optfid, int num)
{
	tN = N_T; xN = 1; yN = 1; zN = N_Z;
	piece_size=N_T; iscomplex=0;
	init_dumpfile(optfid, num);
}

void dump_type_average_spectrum_class :: dump_noflush()
{
	float_type* b1 = (float_type*)NL_SMALLBUFFER1;
	float_type* b2 = b1 + N_T;
	for (size_t cw=0; cw<N_T; cw++)
	{
	 float_type *cb1 = b1+cw;
	 f_complex  *cf  = FIELD+cw;

	 (*cb1)=0;
#ifndef _UNIAXIAL
	 for (size_t cx=0; cx<MY_NX; cx++) for (size_t cy=0; cy<MY_NY; cy++) b1[cw]+=abs2(FIELD[cw+N_T*(cx+N_X*cy)]); 
	 b1[cw] /= N_X*N_Y/(XMAX-XMIN)/(YMAX-YMIN);
#else
	 float_type x1 = HT_PLAN->x_n(MY_NXstart),x2=0; 
	 if (MY_NX>1) x2 = HT_PLAN->x_n(MY_NXstart+1); else {x2=x1; x1=HT_PLAN->x_n(MY_NXstart-1);}
	 float_type f1 = abs2(cf[0]), f2 = 0; 
	 (*cb1) = f1*x1*(x2-x1);
	 for (size_t cx=1; cx<MY_NX; cx++) 
	 {
		 x2 = HT_PLAN->x_n(cx+MY_NXstart); 
		 f1 = abs2(cf[N_T*cx]); 
		 (*cb1)+=f1*x2*(x2-x1);
		 x1 = x2; 
	 }
	 (*cb1)*=XMAX*XMAX/2.0;
#endif
	}
	MPI_Reduce(b1, b2, N_T, MPI_FLOAT_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ISMASTER) dump_write(b2); 
}


dump_type_plasma_axis_class :: dump_type_plasma_axis_class(FILE* optfid, int num)
{
	tN = N_T; xN = 1; yN = 1; zN = N_Z;
	piece_size=N_T; iscomplex=0;
	init_dumpfile(optfid, num);
}


void dump_type_plasma_axis_class :: dump_noflush()
{
#ifdef _UNIAXIAL
	if (!ISMASTER) return;
	int nx = 0; 
	int ny = 0; 
#else
	int ny = N_Y/2 - MY_NYstart; 
	if ((ny<0) || (MY_NY<=ny)) return;
	int nx = N_X/2;
#endif
	float_type* ro = (float_type*)(BIGBUFFER1+N_T);
	fftwt_execute_dft(FFT_BWPLAN_T, (fftwt_complex*)FIELD+N_T*(nx+MY_NX*ny), (fftwt_complex*)BIGBUFFER1); fftwt_Nnormalize(1, BIGBUFFER1);
	calculate_plasmadensity_losses_small(BIGBUFFER1, ro);
	dump_write(ro);
}

void dump_type_Z_axis_class :: dump_noflush()
{
#ifdef _UNIAXIAL
	if (!ISMASTER) return;
	int nx = 0; 
	int ny = 0; 
#else
	int ny = N_Y/2 - MY_NYstart; 
	if ((ny<0) || (MY_NY<=ny)) return;
	int nx = N_X/2;
#endif
	float_type* ro = (float_type*)(BIGBUFFER1+N_T);
	fftwt_execute_dft(FFT_BWPLAN_T, (fftwt_complex*)FIELD+N_T*(nx+MY_NX*ny), (fftwt_complex*)BIGBUFFER1); fftwt_Nnormalize(1, BIGBUFFER1);
	calculate_plasmadensity_losses_small(BIGBUFFER1, ro); for (int nt=0; nt<N_T; nt++) ro[nt]/=NEUTRAL_DENSITY;
	dump_write(ro);
}




dump_type_full_plasma_class :: dump_type_full_plasma_class(FILE* fid, int num)
{
	xN = N_X; yN=N_Y; zN=N_Z; tN=N_T;  
	piece_size = MY_NX*tN*MY_NY;
	iscomplex = 0;
	init_dumpfile(fid, num);
}

void dump_type_full_plasma_class :: dump_noflush()
{
	float_type* b1   = (float_type*)BIGBUFFER1; 
	float_type* nlb1 = (float_type*)NL_SMALLBUFFER1; 
	fftwt_execute(FFT_ALLBWPLAN_T); fftwt_Nnormalize(N_X*MY_NY, BIGBUFFER1);
	//The following code appears somewhat weird, since it takes field values from BIGBUFFER 1
	// and calculates plasma density also in BIGBUFFER1. However, this allows to save the memory space 
	// required and not to use BIGBUFFER2, which contains field x-y Fourier transform and should not be modified. 
	calculate_plasmadensity_losses_small(BIGBUFFER1, nlb1); memcpy(b1, nlb1, N_T*sizeof(float_type)); 
	for (long int n=1; n < MY_NX*MY_NY; n++) calculate_plasmadensity_losses_small(BIGBUFFER1+n*N_T, b1+n*N_T);  
	//calculate_plasmadensity_2float(BIGBUFFER1, b1, MY_NX*MY_NY);
	dump_write(b1); 
}

dump_type_plasma_max_class :: dump_type_plasma_max_class(FILE* optfid, int num)
{
	xN=N_X; yN=N_Y; tN=1; zN=N_Z;
	piece_size = MY_NX*MY_NY;
	iscomplex = 0;
	init_dumpfile(optfid, num);
}

void dump_type_plasma_max_class :: dump_noflush()
{
	float_type* b1  = (float_type*)BIGBUFFER1;
	float_type* nlb = (float_type*)NL_SMALLBUFFER1; 

	fftwt_execute(FFT_ALLBWPLAN_T); fftwt_Nnormalize(MY_NX*MY_NY, BIGBUFFER1);

	calculate_maxplasmadensity_2float(BIGBUFFER1, b1, ((size_t)(MY_NX*MY_NY)));
	
	
	/*for (int ny=0; ny<MY_NY; ny++)
	for (int nx=0; nx<N_X;   nx++)
	{
		float_type maxro = 0;
		int ofs1 = (nx+N_X*ny);
		calculate_plasmadensity_losses_small(BIGBUFFER1 + N_T*ofs1, nlb); 
		for (int nt=0; nt<N_T; nt++) maxro = max(maxro, nlb[nt]);
		b1[ofs1] = maxro;
	}*/
	dump_write(b1);
}


dump_type_ysection_plasma_class :: dump_type_ysection_plasma_class(FILE* optfid, int num)
{
	tN=N_T; xN = N_X; yN=1; zN=N_Z;
	piece_size=tN*MY_NX; iscomplex=0;
	init_dumpfile(optfid, num);
	ny = load_namednumint(optfid, "DUMP_NY", num, true, N_Y/2);
	ny -= MY_NYstart; if (0 <= ny && ny < MY_NY) mydump=true; else mydump=false;
}


void dump_type_ysection_plasma_class :: dump_noflush()
{
	if (mydump) 
	{
		float_type* b   = (float_type*)BIGBUFFER1;
		float_type* nlb = (float_type*)NL_SMALLBUFFER1; 
		fftwt_execute(fft_bwplan_sectiont); fftwt_Nnormalize(MY_NX, BIGBUFFER1);
		//calculate_plasmadensity_losses_small(BIGBUFFER1, nlb); memcpy(b, nlb, N_T*sizeof(float_type)); 
		//for (int nx=1; nx<N_X; nx++) calculate_plasmadensity_losses_small(BIGBUFFER1+N_T*nx, b+N_T*nx);
		calculate_plasmadensity_2float(BIGBUFFER1, b, MY_NX);
		dump_write(b);
	}
}

dump_type_ysection_plasma_max_class :: dump_type_ysection_plasma_max_class(FILE* optfid, int num)
{
	tN=1; xN = N_X; yN=1; zN=N_Z;
	piece_size=tN*MY_NX; iscomplex=0;
	init_dumpfile(optfid, num);
	init_ny(optfid, num);
	fft_bwplan_sectiont = fftwt_plan_many_dft(1, &N_T, N_X, (fftwt_complex*)FIELD+N_T*MY_NX*ny, NULL, 1, N_T,  (fftwt_complex*)BIGBUFFER1, NULL, 1, N_T, FFTW_BACKWARD, FFTW_FLAG); 


}


void dump_type_ysection_plasma_max_class :: dump_noflush()
{
	if (mydump) 
	{
		float_type* b   = (float_type*)BIGBUFFER1;
		float_type* nlb = (float_type*)NL_SMALLBUFFER1;  
		fftwt_execute(fft_bwplan_sectiont); fftwt_Nnormalize(MY_NX, BIGBUFFER1);
		for (int nx=0; nx<N_X; nx++) 
		{
			float_type maxro = 0;
			calculate_plasmadensity_losses_small(BIGBUFFER1+N_T*nx, nlb);
			for (int nt=0; nt<N_T; nt++) maxro = max(maxro, nlb[nt]);
			b[nx]=maxro;
		}
		dump_write(b);
	}
}

