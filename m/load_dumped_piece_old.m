function [P, Z, Nt, tmin, tmax, xnet, ynet,nodes_x, nodes_y] = load_dumped_piece_old(filename)

fid = fopen(filename,'rb');
Z  = fread(fid, 1, 'double');
Nt = fread(fid, 1, 'int');
tmin = fread(fid, 1, 'double');
tmax = fread(fid, 1, 'double');

Nx = fread(fid, 1, 'int');
xnet = fread(fid, Nx, 'double');
Ny = fread(fid, 1, 'int');
ynet = fread(fid, Ny, 'double');

nodes_x = fread(fid, 2, 'int');
nodes_y = fread(fid, 2, 'int');

myrankx = fread(fid, 1, 'int');
myranky = fread(fid, 1, 'int'); 

myNx = nodes_x(2)-nodes_x(1)-1;
myNy = nodes_y(2)-nodes_y(1)-1;

P = fread(fid, 2*Nt*myNy*myNy, 'double');
P = reshape(complex(P(1:2:end), P(2:2:end)), Nt, myNx, myNy);


fclose(fid);

