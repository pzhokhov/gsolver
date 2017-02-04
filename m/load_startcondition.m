function [A, tmin,tmax, xmin,xmax,ymin,ymax, znet, Vg, omega0, omega, wavenum0, wavenum] = load_startcondition(filename)
fname = ls(filename)
fid = fopen(fname,'rb')

omega0 = fread(fid, 1, 'double');
wavenum0 = fread(fid, 2, 'double'); wavenum0 = complex(wavenum0(1), wavenum0(2));
Vg       = fread(fid, 1, 'double'); 

Nt = fread(fid, 1, 'int');
tmin = fread(fid, 1, 'double');
tmax = fread(fid, 1, 'double');

omega  = fread(fid, Nt, 'double');
wavenum     = fread(fid, 2*Nt, 'double');

wavenum    = complex(   wavenum(1:2:end),    wavenum(2:2:end));

Nx   = fread(fid, 1, 'int');
xmin = fread(fid, 1, 'double');
xmax = fread(fid, 1, 'double');

Ny   = fread(fid, 1, 'int');
ymin = fread(fid, 1, 'double');
ymax = fread(fid, 1, 'double');

Nz   = fread(fid, 1, 'int');
znet = fread(fid, Nz, 'double');

A0 = fread(fid, 2*Nt*Nx*Ny, 'double');
A  = complex(A0(1:2:end), A0(2:2:end));
clear A0;
A = reshape(A, Nt, Nx, Ny);

fclose(fid);
