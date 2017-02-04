function [E,w,phix,phiy]= xyangle_interpolation_file(filename, phix, phiy)
% function [E,w,phix]=xangle_interpolation_file(filename, nz)
%
% This function loads field from gsolver dump file named filename at z
% point number nz and runs xangle_interpolation for it.
%
% See also xangle_interpolation

%[A, Nt, tmin, tmax, Nx, xmin, xmax, Ny, ymin, ymax, znet, omega0, Vg, omega, wavenum0, wavenum] = load_startsectiony(filename,ny);
[A, tmin, tmax, xmin, xmax,  ymin, ymax, znet, Vg, omega0, omega, wavenum0, wavenum] = load_startcondition(filename);
if (nargin < 2)
	phix = -0.5 : 0.005 : 0.5;
end;
if (nargin < 3)
	phiy = -0.2 : 0.002 : 0.2;
end;

[E, w, phix, phiy] = xyangle_interpolation(A,omega,wavenum, xmin, xmax, ymin, ymax, phix, phiy);
