function [P, Z, Nt, tmin, tmax, xnet, ynet,myNy,myNystart] = load_dumped_piece(filename)

fid = fopen(filename,'rb');
Z  = fread(fid, 1, 'double');
Nt = fread(fid, 1, 'int');
tmin = fread(fid, 1, 'double');
tmax = fread(fid, 1, 'double');

Nx = fread(fid, 1, 'int');
xnet = fread(fid, Nx, 'double');
Ny = fread(fid, 1, 'int');
ynet = fread(fid, Ny, 'double');

myNy = fread(fid, 1, 'int');
myNystart = fread(fid, 1, 'int');

P = fread(fid, 2*Nt*Nx*myNy, 'double');
P = reshape(complex(P(1:2:end), P(2:2:end)), Nt, Nx, myNy);


fclose(fid);

