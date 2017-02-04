function [A, Zout, Nt, tmin, tmax,xnet, ynet] = load_dump_sectiony(filemask, ny)
files = dir(filemask);
if isempty(files)
    error('No files match requred name mask.');
end;
[P, Z, Nt, tmin, tmax, xnet, ynet, myNy] = load_dumped_piece(files(1).name);

process_Ny = length(ynet)./(myNy);
N  = length(files)/process_Ny;
if N ~= floor(N)
    error('Total number of files is not proportional to the number of files in one step. Possibly wrong filemask');
end;

A    = complex(zeros(Nt, length(xnet), N));
if (nargin < 2); ny=length(ynet)/2; end;
Zout = zeros(N,1); Zout(1) = Z;

for k = 1 : N
for n = 1 : process_Ny;
    [P, Zout(k), Nt, tmin, tmax, xnet, ynet, myNy, myNystart] = load_dumped_piece(files(n+(k-1)*process_Ny).name);
    if (ny < myNystart+1 || myNystart+myNy > ny); continue; end;
    for nx = 1 : length(xnet)
    for nt = 1 : Nt;
        A(nt, nx, k) = P(nt, nx, ny-myNystart);
    end;
    end;
    disp(sprintf('File %s loaded', files(n+(k-1)*process_Ny).name));    
end;
end;