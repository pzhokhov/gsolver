function [Imax, z] = getGDmaxI(filemask, step)
%getGDmaxI
%
% This function extracts maximum intensity from GNLSsolver dumps
% call format:
%
% [Imax,z] = getGDmaxI(filemask, step)
%
% filemask - mask of dump file names
% step     - this function takes every 'step' dump in z direction
%            assumed to be 1 if ommitted

files = dir(filemask);

if (nargin < 1)
    error('Too few arguements specified!');
elseif (nargin == 1)
    step = 1;
end;

[P, z0, Nt, tmin, tmax, xnet, ynet, nodes_x, nodes_y] = load_dumped_piece(files(1).name);

process_Nx = length(xnet)./(nodes_x(2)-nodes_x(1)-1);
process_Ny = length(ynet)./(nodes_y(2)-nodes_y(1)-1);
N  = length(files)/process_Nx/process_Ny;
if N ~= floor(N)
    error('Total number of files is not proportional to the number of files in one step. Possibly wrong filemask');
end;

z    =  zeros(N,1);
Imax =  zeros(N,1);


for i = 1 : step*length(files)/N : size(files,1);
    name = [files(i).name(1:end-10), '*'] ; [A,z(i)] = load_dump(name);
    Imax(i) = max(abs(A(:))).^2;
end;