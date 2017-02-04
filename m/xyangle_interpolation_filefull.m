function [E,w,phix,phiy]= xyangle_interpolation_filefull(filename)
% function xyangle_interpolation_filefull(filename, nz)
%
% This function loads field from gsolver final dump file named filename and
% runs xyangle_interpolation for it.
%
% See also xyangle_interpolation

%[A, Nt, tmin, tmax, Nx, xmin, xmax, Ny, ymin, ymax, znet, omega0, Vg, omega, wavenum0, wavenum] = load_startsectiony(filename,ny);
[A,Nt, tmin,tmax, Nx, xmin,xmax, Ny,ymin,ymax, znet, omega0, omega, wavenum0, wavenum, Vg] = load_startcondition(filename);
[E,w,phix,phiy] = xyangle_interpolation(A,omega,wavenum, xmin, xmax, ymin, ymax);