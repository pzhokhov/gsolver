Ok so this is Sasha's README file - taking notes about how to use peter's filamentation code. 

Waveion is the fdtd code for ionization off a surface - includes Peter's new and fancy ionization equations (i.e. non-Keldysh to allow for sub-cycle pulse changes).

The main filamentation code is contained in gsolver and the matlab code to load the results is in m. Bin contains the output of runs and run configurations (pulse power, output results, etc.)

Physical processes taken into account: 
1) single mode Raman (with given frequency and decay time). 
2) instantaneous Kerr; self-steepening
3) Plasma; Keldysh ionization, decay with constant time, avalanche ionization (goverend by collision time), plasma absorption and defocusing.

All of the code needs gsl, gslcblas, cuda compiler, openmpi, and fftw3 with openmpi flags to run.


