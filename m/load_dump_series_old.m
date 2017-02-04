function [A, Zout, Nt, tmin, tmax,xnet, ynet] = load_dump_series_old(filemask)

files = dir(filemask);
if isempty(files)
    error('No files match requred name mask.');
end;
[P, Z, Nt, tmin, tmax, xnet, ynet, nodes_x, nodes_y] = load_dumped_piece_old(files(1).name);

process_Nx = length(xnet)./(nodes_x(2)-nodes_x(1)-1);
process_Ny = length(ynet)./(nodes_y(2)-nodes_y(1)-1);
N  = length(files)/process_Nx/process_Ny;
if N ~= floor(N)
    error('Total number of files is not proportional to the number of files in one step. Possibly wrong filemask');
end;

A    = complex(zeros(Nt, length(xnet), length(ynet), N));
Zout = zeros(N,1); Zout(1) = Z;

for k = 1 : N
for n = 1 : process_Nx*process_Ny;
    [P, Zout(k), Nt, tmin, tmax, xnet, ynet, nodes_x, nodes_y] = load_dumped_piece_old(files(n+(k-1)*process_Nx*process_Ny).name);
    for ny = nodes_y(1)+2 : nodes_y(2)
    if ny == 0; continue; end;
    for nx = nodes_x(1)+2 : nodes_x(2)
    if nx == 0; continue; end;
    for nt = 1 : Nt;
        A(nt, nx, ny, k) = P(nt, nx-nodes_x(1)-1, ny-nodes_y(1)-1);
    end;
    end;
    end;
    disp(sprintf('File %s loaded', files(n+(k-1)*process_Nx*process_Ny).name));    
end;
end;