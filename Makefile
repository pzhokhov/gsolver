ofiles = solver_main.o process_input.o propagator.o ionization.o dumps.o dump_types.o extmath.o \
 	  pulse_gg.o pulse_beamsgg.o pulse_beamsprobegg.o pulse_beamsprobeextgg.o pulse_beamsprobevgg.o pulse_beams_delayedprobegg.o        \
	  pulse_beams_delayedprobegg_randphase.o pulse_customspectrumg.o fhatha.o fhatha_cuda.o

output_filename = gsolver_gas

hfiles  = solver.h extmath.h dumps.h pulseshapes.h ionization.h SIconstants.h
cuda_hfiles = cuda_extmath.h


CC	     = mpicxx
NVCC	     = nvcc -I /usr/include/openmpi -arch=sm_20
#NVCC	     = nvcc -I /usr/include/openmpi -arch=sm_20
#nvcc -arch=sm_20 --ptxas-options=-v -I${MPI_INCLUDE}
output_dir   = bin/


gsolver_gas : CPPFLAGS = -O3 -g -DPPT_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -DTRANSVERSE_DIMENSIONS=1 -DPLASMA_DISPERSION -DTHIRD_HARMONICS -DMULTI_LEVEL_IONIZATION -DFFTW_FLAG=FFTW_MEASURE 
gsolver_gas_qdht : CPPFLAGS = -O3 -g -DPPT_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -D_UNIAXIAL -DPLASMA_DISPERSION -DTHIRD_HARMONICS -DMULTI_LEVEL_IONIZATION -DFFTW_FLAG=FFTW_MEASURE -DQDHT -DTRANSVERSE_DIMENSIONS=1


gsolver_gas_qdht_cuda : CPPFLAGS = -O0 -g -DPPT_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -D_UNIAXIAL -DPLASMA_DISPERSION -DTHIRD_HARMONICS -DMULTI_LEVEL_IONIZATION -DFFTW_FLAG=FFTW_MEASURE -DQDHT -DTRANSVERSE_DIMENSIONS=1 -D_ENABLE_CUDA

# gsolver_solid_qdht is 1D axially symmetric beam using quasi-discrete Henkel transfsorm; main diff between this and ux is DGDHT flag at the end.LAST FLAG ADDED POST SASHA

gsolver_solid_qdht : CPPFLAGS = -O3 -g -DKELDYSH_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -D_UNIAXIAL -DPLASMA_DISPERSION  -DFFTW_FLAG=FFTW_MEASURE -DQDHT -DTRANSVERSE_DIMENSIONS=1 

 
gsolver_solid_ux : CPPFLAGS = -O3 -g -DKELDYSH_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -DPLASMA_DISPERSION -DFFTW_FLAG=FFTW_MEASURE -DTRANSVERSE_DIMENSIONS=1
gsolver_solid_1d : CPPFLAGS = -O3 -g -DKELDYSH_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -DPLASMA_DISPERSION -DFFTW_FLAG=FFTW_MEASURE -DNO_DIFFRACTION
gsolver_solid : CPPFLAGS = -O3 -g -DKELDYSH_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -DPLASMA_DISPERSION -DFFTW_FLAG=FFTW_MEASURE -DTRANSVERSE_DIMENSIONS=2
gsolver_solid_qdht : CPPFLAGS = -O3 -g -DKELDYSH_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -DPLASMA_DISPERSION -DFFTW_FLAG=FFTW_ESTIMATE -DTRANSVERSE_DIMENSIONS=1 -DQDHT -D_UNIAXIAL

gsolver_bloch : CPPFLAGS = -O0 -g -DPPT_IONIZATION -DNO_PLASMA_RESPONSE -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -DFFTW_FLAG=FFTW_MEASURE -DTRANSVERSE_DIMENSIONS=0 -DBLOCH_RESPONSE


gsolver_SS : CPPFLAGS = -O3 -g -DPPT_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -D_UNIAXIAL -DQDHT -DNO_PLASMARESPONSE -DNO_SPACE_TIME_FOCSING -DFFTW_FLAG=FFTW_MEASURE 

gsolver_SS_STF : CPPFLAGS = -O0 -g -DPPT_IONIZATION -D_EXTMATH_DOUBLE -D_FILE_OFFSET_BITS=64 -D_ODE_RK4 -D_UNIAXIAL -DQDHT -DNO_PLASMARESPONS -DTHIRD_HARMONICS -DFFTW_FLAG=FFTW_MEASURE 

NVCCFLAGS           = $(CPPFLAGS)


${MAKECMDGOALS} : $(ofiles)
	$(CC) $(ofiles) $(CPPFLAGS) -lfftw3_mpi -lfftw3 -lgsl -lgslcblas -o $(output_dir)/$@
	cp $(output_dir)/$@ .
	
%.o : %.cpp $(hfiles)
	$(CC) -c $< -o $@ $(NVCCFLAGS)

%.o : %.cu $(hfiles) $(cuda_hfiles)
	$(CC) -x c++ -c $< -o $@ $(NVCCFLAGS)



.DEFAULT_GOAL=gsolver_gas_qdht
.PHONY : clean
clean :
	rm -f $(ofiles)

