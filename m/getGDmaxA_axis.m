function [Imax, z] = getGDmaxI_axis(filemask, step)
%getGDmaxI_axis
%
% This function extracts maximum on-axis intensity from GNLSsolver dumps
% call format:
%
% [Imax,z] = getGDmaxI_axis(filemask, step)
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

istep = step*length(files)/N;
z    =  zeros(floor(N/istep),1);
Imax =  zeros(floor(N/istep),1);


hwait = waitbar(0,'Loading data, please wait...');

for i = 1 : istep : size(files,1)
    name = [files(i).name(1:end-10), '*'] ; [A,z((i-1)/istep+1)] = load_dump(name);
    Imax((i-1)/istep+1) = max(abs(A(:, size(A,2)/2, size(A,3)/2)));
    waitbar(i/size(files,1), hwait);
end;

close(hwait);