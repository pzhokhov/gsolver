This README is to guide you through what the config file (i.e. fil_<Material>.txt) does, since the config file itself does not support comments. 

INITIAL PARAMETERS

lambda_V = central wavelength of the pulse in vacuum

N_T = number of poitns in time

T_MIN/T_MAX = temporal window from T_MIN to T_MAX (assuming pulse is at T = 0). each time step = (T_MAX - T_MIN)/N_T

DISPERSION_FILE = file from which refractive index vs lambda is read - see write_refractive_index in the m directory.

N_X, N_Y, N_R = # of poitns in transverse directions

X_MIN/MAX, ... = size of spatial window

Z_MIN/MAX = propogation window

ZNET_TYPE = couple of options here: e = equidistant, each step is the same size. p3 = polynomial of the 3rd degree, there are more points in the center than the edge. 

PULSE_SHAPE = typically gg, gaussian in space and time. 

DURATION_FWHM = self explanatory, duration at the full-width at half maximum.

DIAMTER_FWHM = self-explanatory...

WAVEFRONT_RADIUS = focusing distance, focus with respect to Z_MIN. In SI units, make this bigger if the code is running too slow. 

ENERGY = In J

NOISE_LEVEL = Random white noise level

CHIRP_SPECTRAL = SI units

NONLIN_REFRINDEX = m^2/W (SI units)

RAMAN_FRACTION = fraction of delayed response (1 = full delay, 0 = instantaneous nonlinear response)

TAU_RAMAN/OMEGA_RAMAN: Raman equation, get back to this.

NEUTRAL_DENSITY = Density of atoms

RECOMBINATION_TAU = recombination time, trapping of free electrons. Leave as is generally.

AMBIENT_CARRIERS = ambient free electrons, usually 0.

COLLISION_TAU = single electron time between collisions. Can look this up.

IONIZATION_POTENTIAL = taken to be same as bandgap - can be looked up.


DUMPS

DUMPS_N = number of different outputs, here are the different types:

"fa": The abs value of the "fa" dump squared gives the power in units of 10^(12) W/cm^2. 

"f": Flux - This is equal to the integral of I over time.

"pm": plasma max - maximum density of plasma with respect to time.

"full": this is the field in each pt of the beam at each point in time, every Nz steps saves in one big file. 

"mI": Maximum Intensity, needs no further explanation. 
