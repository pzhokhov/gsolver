function [E,w,phix]= xangle_interpolation_file(filename, nz)
% function [E,w,phix]=xangle_interpolation_file(filename, nz)
%
% This function loads field from gsolver dump file named filename at z
% point number nz and runs xangle_interpolation for it.
%
% See also xangle_interpolation

%[A, Nt, tmin, tmax, Nx, xmin, xmax, Ny, ymin, ymax, znet, omega0, Vg, omega, wavenum0, wavenum] = load_startsectiony(filename,ny);
   [A, Nt, tmin, tmax, Nx, xmin, xmax,                 znet, Vg, omega0, wavenum0, omega, wavenum] = load_sectiondump(filename, nz);
[E, w, phix] = xangle_interpolation(A,omega,wavenum, xmin, xmax);
