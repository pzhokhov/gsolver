function [A,Nt, tmin,tmax, xnet, ynet, znet, omega0, Vg, omega, wavenum0, wavenum] = load_startsectiony(filename,ny)

fname = ls(filename);
fid = fopen(fname,'rb');

omega0 = fread(fid, 1, 'float');
wavenum0 = fread(fid, 2, 'float'); wavenum0 = complex(wavenum0(1), wavenum0(2));
Vg       = fread(fid, 1, 'float'); 

Nt = fread(fid, 1, 'int');
tmin = fread(fid, 1, 'float');
tmax = fread(fid, 1, 'float');

omega  = fread(fid, Nt, 'float');
wavenum     = fread(fid, 2*Nt, 'float');

wavenum    = complex(   wavenum(1:2:end),    wavenum(2:2:end));

Nx   = fread(fid, 1, 'int');
xmin = fread(fid, 1, 'float');
xmax = fread(fid, 1, 'float');

Ny   = fread(fid, 1, 'int');
ymin = fread(fid, 1, 'float');
ymax = fread(fid, 1, 'float');

Nz   = fread(fid, 1, 'int');
znet = fread(fid, Nz, 'float');

fseek(fid, 2*(Nt*Nx*ny+1)*4, 'cof');

A0 = fread(fid, 2*Nt*Nx, 'float');
A = complex(A0(1:2:end), A0(2:2:end));
A = reshape(A, Nt, Nx);

fclose(fid);
