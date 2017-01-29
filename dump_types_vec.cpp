#include "dumps.h"

#include "solver.h"

dump_type_empty_class :: ~dump_type_empty_class()
{
	fclose(dumpfid);
}


void dump_type_empty_class :: init_dumpfile(FILE* optfid, int i)
{
	char filenamebuf[200];
	load_namednumstringn(optfid, "DUMP_FILENAME", filenamebuf, i, 200, false);
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
#ifdef __ENABLE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
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
	size_t ofs = header_ofs; 
	if (iscomplex) ofs += fileofs()*sizeof(f_complex); else ofs += fileofs()*sizeof(float_type);
	fseeko(dumpfid, ofs, SEEK_SET);
	if (iscomplex) fwrite(buf, sizeof(f_complex) , piece_size, dumpfid);
	else           fwrite(buf, sizeof(float_type), piece_size, dumpfid);
}




dump_type_full_class :: dump_type_full_class(FILE* fid, int num)
{
	xN = N_X; yN=N_Y; zN=N_Z; tN=2*N_T;  
	piece_size = xN*tN*MY_NY;
	iscomplex = 1; 
	init_dumpfile(fid, num);
}


size_t dump_type_full_class :: fileofs()
{
	return tN*xN*(yN*((size_t)n_Z) + ((size_t)MY_NYstart));
}

void dump_type_full_class :: dump_noflush()
{
	
	fftw(FFT_BWPLAN_T, N_X*MY_NY, (fftwt_complex*)FIELD,     2,2*N_T, (fftwt_complex*)BIGBUFFER1,     2,2*N_T); 
	fftw(FFT_BWPLAN_T, N_X*MY_NY, (fftwt_complex*)(FIELD+1), 2,2*N_T, (fftwt_complex*)(BIGBUFFER1+1), 2,2*N_T); fftwt_Nnormalize(2*N_X*MY_NY, BIGBUFFER1);
	dump_write(BIGBUFFER1); 
}


dump_type_maxI_class :: dump_type_maxI_class(FILE* fid, int num)
{
	xN = N_X; yN=N_Y; zN=N_Z; tN=1; 
	piece_size = xN*MY_NY;
	iscomplex = 0; 
	init_dumpfile(fid, num);

	wfilter = (f_complex*)malloc_ch(2*sizeof(f_complex)*N_T); 
	char filtertype[100];  load_namednumstringn(fid, "FILTER_TYPE", filtertype, num, 100, "true", "NO");
	create_wfilter(wfilter, filtertype);
}

dump_type_maxI_class :: ~dump_type_maxI_class() 
{
	free(wfilter);
}

size_t dump_type_maxI_class :: fileofs()
{
	return xN*(yN*n_Z+MY_NYstart);
}

void dump_type_maxI_class :: dump_noflush()
{
	float_type* buf = (float_type*)BIGBUFFER1;
	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<N_X; nx++)
	{
		int ofs = 2*N_T*(nx+N_X*ny);
		for(int nw=0; nw<2*N_T; nw++) NL_SMALLBUFFER1[nw]=FIELD[ofs+nw]*wfilter[nw];
		fftw(FFT_BWPLAN_T, 2, (fftwt_complex*)NL_SMALLBUFFER1, 2, 1, (fftwt_complex*)(BIGBUFFER1+ofs), 2, 1); 
	}
    fftwt_Nnormalize(2*N_X*MY_NY, BIGBUFFER1);

	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<N_X; nx++)
	{
		float_type maxI=0;
		for(int nt=0; nt<N_T; nt++) 
		{
			int ofs = 2*(nt + N_T*(nx+N_X*ny)); maxI = max(maxI, abs2(BIGBUFFER1[ofs])+abs2(BIGBUFFER1[ofs+1]));
		}
		buf[nx+N_X*ny] = maxI;
	}
	dump_write(buf);
}

void dump_type_flux_class :: dump_noflush()
{
	// TODO: Compare results for frequency- and time-domain integration
	float_type* buf = (float_type*)BIGBUFFER1;
	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<N_X; nx++)
	{
		int ofs = 2*N_T*(nx+N_X*ny);
		for(int nw=0; nw<2*N_T; nw++) NL_SMALLBUFFER1[nw]=FIELD[ofs+nw]*wfilter[nw];
		fftw(FFT_BWPLAN_T, 1, (fftwt_complex*)NL_SMALLBUFFER1, 1, N_T, (fftwt_complex*)(BIGBUFFER1+ofs), 1, N_T); 
	}
    fftwt_Nnormalize(N_X*MY_NY, BIGBUFFER1);

	for(int ny=0; ny<MY_NY; ny++) for(int nx=0; nx<N_X; nx++)
	{
		float_type f=0;
		for(int nt=0; nt<2*N_T; nt++) 
		{
			int ofs = nt + N_T*(nx+N_X*ny); f+=abs2(BIGBUFFER1[ofs]);
		}
		buf[nx+N_X*ny] = f*TSTEP;
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
	tN=2*N_T; xN = N_X; yN=1; zN=N_Z;
	piece_size=tN*xN; iscomplex=1;
	init_dumpfile(optfid, num);
	init_ny(optfid, num);
}

size_t dump_type_ysection_class :: fileofs()
{
	return piece_size*n_Z;
}

void dump_type_ysection_class :: dump_noflush()
{
	if (mydump) 
	{
		fftw(FFT_BWPLAN_T, N_X, (fftwt_complex*)FIELD+2*N_T*N_X*ny,   2,2*N_T, (fftwt_complex*)BIGBUFFER1,   2,2*N_T);
		fftw(FFT_BWPLAN_T, N_X, (fftwt_complex*)FIELD+2*N_T*N_X*ny+1, 2,2*N_T, (fftwt_complex*)BIGBUFFER1+1, 2,2*N_T);

		fftwt_Nnormalize(N_X, BIGBUFFER1);
		dump_write(BIGBUFFER1);
	}
}


dump_type_ysection_maxI_class :: dump_type_ysection_maxI_class(FILE* optfid, int num)
{
	throw "TODO: vector version";
	tN=1; xN = N_X; yN=1; zN=N_Z;
	piece_size=xN; iscomplex=0;
	init_dumpfile(optfid, num);
	init_ny(optfid, num);

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
		fftw(FFT_BWPLAN_T, N_X, (fftwt_complex*)FIELD+N_T*N_X*ny, 1,N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(N_X, BIGBUFFER1);
		for (int nx=0; nx<N_X; nx++) 
		{
			float_type maxI = 0;
			for (int nt=0; nt<N_T; nt++) maxI=max(maxI, abs2(BIGBUFFER1[nt+N_T*nx]));
			buf[nx] = maxI;
		}
		dump_write(buf);
	}
}

void dump_type_ysection_flux_class :: dump_noflush()
{
	if (mydump) 
	{
		float_type* buf = (float_type*)BIGBUFFER1;
		fftw(FFT_BWPLAN_T, N_X, (fftwt_complex*)FIELD+N_T*N_X*ny, 1,N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(N_X, BIGBUFFER1);
		for (int nx=0; nx<N_X; nx++) 
		{
			float_type f = 0;
			for (int nt=0; nt<N_T; nt++) f+=abs2(BIGBUFFER1[nt+N_T*nx]);
			buf[nx] = f*TSTEP;
		}
		dump_write(buf);
	}
}



dump_type_full_plasma_class :: dump_type_full_plasma_class(FILE* fid, int num)
{
	throw "TODO: vector version";
	xN = N_X; yN=N_Y; zN=N_Z; tN=N_T;  
	piece_size = xN*tN*MY_NY;
	iscomplex = 0;
	init_dumpfile(fid, num);
}

void dump_type_full_plasma_class :: dump_noflush()
{
	float_type* b1   = (float_type*)BIGBUFFER1; 
	float_type* nlb1 = (float_type*)NL_SMALLBUFFER1; 
	fftw(FFT_BWPLAN_T, N_X*MY_NY, (fftwt_complex*)FIELD, 1, N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(N_X*MY_NY, BIGBUFFER1);
	//The following code appears somewhat weird, since it takes field values from BIGBUFFER 1
	// and calculates plasma density also in BIGBUFFER1. However, this allows to save the memory space 
	// required and not to use BIGBUFFER2, which contains field x-y Fourier transform and should not be modified. 
	calculate_plasmadensity_small_2float(BIGBUFFER1, nlb1); memcpy(b1, nlb1, N_T*sizeof(float_type)); 
	calculate_plasmadensity_2float(BIGBUFFER1+N_T, b1+N_T, N_X*MY_NY-1);  
	dump_write(BIGBUFFER2); 
}

dump_type_plasma_max_class :: dump_type_plasma_max_class(FILE* optfid, int num)
{
	throw "TODO: vector version";
	xN=N_X; yN=N_Y; tN=1; zN=N_Z;
	piece_size = xN*tN*MY_NY;
	iscomplex = 0;
	init_dumpfile(optfid, num);
}

void dump_type_plasma_max_class :: dump_noflush()
{
	throw "TODO: vector version";
	float_type* b1  = (float_type*)BIGBUFFER1;
	float_type* nlb = (float_type*)NL_SMALLBUFFER1; 

	fftw(FFT_BWPLAN_T, N_X*MY_NY, (fftwt_complex*)FIELD, 1, N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(N_X*MY_NY, BIGBUFFER1);
	
	calculate_maxplasmadensity_2float(BIGBUFFER1, (float_type*)BIGBUFFER1, N_X*MY_NY);
	/*for (int ny=0; ny<MY_NY; ny++)
	for (int nx=0; nx<N_X;   nx++)
	{
		float_type maxro = 0;
		int ofs1 = (nx+N_X*ny);
		calculate_plasmadensity_small_2float(BIGBUFFER1 + N_T*ofs1, nlb); 
		for (int nt=0; nt<N_T; nt++) maxro = max(maxro, nlb[nt]);
		b1[ofs1] = maxro;
	}*/
	dump_write(b1);
}

dump_type_ysection_plasma_class :: dump_type_ysection_plasma_class(FILE* optfid, int num)
{
	throw "TODO: vector version";
	tN=N_T; xN = N_X; yN=1; zN=N_Z;
	piece_size=tN*xN; iscomplex=0;
	init_dumpfile(optfid, num);
	ny = load_namednumint(optfid, "DUMP_NY", num, true, N_Y/2);
	ny -= MY_NYstart; if (0 <= ny && ny < MY_NY) mydump=true; else mydump=false;
}


void dump_type_ysection_plasma_class :: dump_noflush()
{

	if (mydump) 
	{
		float_type* b   = (float_type*)BIGBUFFER1;
		//float_type* nlb = (float_type*)NL_SMALLBUFFER1; 
		fftw(FFT_BWPLAN_T, N_X, (fftwt_complex*)FIELD+N_T*N_X*ny, 1,N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(N_X, BIGBUFFER1);
		//calculate_plasmadensity_small_2float(BIGBUFFER1, nlb); memcpy(b, nlb, N_T*sizeof(float_type)); 
		calculate_plasmadensity_2float(BIGBUFFER1, b, N_X);
		dump_write(b);
	}
}

dump_type_ysection_plasma_max_class :: dump_type_ysection_plasma_max_class(FILE* optfid, int num)
{
	throw "TODO: vector version";
	tN=1; xN = N_X; yN=1; zN=N_Z;
	piece_size=tN*xN; iscomplex=0;
	init_dumpfile(optfid, num);
	init_ny(optfid, num);
}


void dump_type_ysection_plasma_max_class :: dump_noflush()
{
	if (mydump) 
	{
		float_type* b   = (float_type*)BIGBUFFER1;
		float_type* nlb = (float_type*)NL_SMALLBUFFER1;  
		fftw(FFT_BWPLAN_T, N_X, (fftwt_complex*)FIELD+N_T*N_X*ny, 1,N_T, (fftwt_complex*)BIGBUFFER1, 1,N_T); fftwt_Nnormalize(N_X, BIGBUFFER1);
		for (int nx=0; nx<N_X; nx++) 
		{
			float_type maxro = 0;
			calculate_plasmadensity_small_2float(BIGBUFFER1+N_T*nx, nlb);
			for (int nt=0; nt<N_T; nt++) maxro = max(maxro, nlb[nt]);
			b[nx]=maxro;
		}
		dump_write(b);
	}
}
