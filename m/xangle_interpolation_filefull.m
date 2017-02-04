function [E,phix]= xangle_interpolation_filefull(filename)
% function xangle_interpolation_file(filename, nz)
%
% This function loads field from gsolver final dump file named filename
% runs xangle_interpolation for it.
%
% See also xangle_interpolation

%[A, Nt, tmin, tmax, Nx, xmin, xmax, Ny, ymin, ymax, znet, omega0, Vg, omega, wavenum0, wavenum] = load_startsectiony(filename,ny);
[A,Nt, tmin,tmax, Nx, xmin,xmax, Ny,ymin,ymax, znet, omega0, omega, wavenum0, wavenum, Vg] = load_startcondition(filename);
[E, phix] = xangle_interpolation(A,omega,wavenum, xmin, xmax);