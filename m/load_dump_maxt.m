function [A, Z, Nt, tmin, tmax,xnet, ynet] = load_dumpmax(filemask)
A = [];
files = dir(filemask);
if isempty(files)
    error('No files match requred name mask.');
end;
[P, Z, Nt, tmin, tmax, xnet, ynet, nodes_x, nodes_y] = load_dumped_piece(files(1).name);
A = complex(zeros(length(xnet), length(ynet)));

hwait = [];
if ((nodes_y(2)-nodes_y(1))*(nodes_x(2)-nodes_x(1))*Nt > 1e8)
    hwait = waitbar(0, sprintf('Loading file %s, please wait', files(1).name));
end;
    for ny = nodes_y(1)+1 : nodes_y(2)
    if ny == 0 continue; end;
    for nx = nodes_x(1)+1 : nodes_x(2)
    if nx == 0 continue; end;
    for nt = 1 : Nt;
        A(nx, ny) = max(A(nx,ny), abs(P(nt+(nx-nodes_x(1)-1)*Nt+(ny-nodes_y(1)-1)*Nt*(nodes_x(2)-nodes_x(1)+1))));
    end;
    if (~isempty(hwait))
        waitbar((nx + ny*(nodes_x(2)-nodes_x(1)-1))/((nodes_x(2)-nodes_x(1)-1)*(nodes_y(2)-nodes_y(1)-1)), hwait);
    end;
    end;
    end;
if ~isempty(hwait)
    close(hwait);
end;
disp(sprintf('File %s loaded', files(1).name));    
for n = 2 : length(files);
    [P, Z, Nt, tmin, tmax, xnet, ynet, nodes_x, nodes_y] = load_dumped_piece(files(n).name);
    for ny = nodes_y(1)+1 : nodes_y(2)
    if ny == 0 continue; end;
    for nx = nodes_x(1)+1 : nodes_x(2)
    if nx == 0 continue; end;
    for nt = 1 : Nt;
        A(nx, ny) = max(A(nx,ny), abs(P(nt, nx-nodes_x(1), ny-nodes_y(1))));
    end;
    if (~isempty(hwait))
        waitbar((nx + ny*(nodes_x(2)-nodes_x(1)-1))/((nodes_x(2)-nodes_x(1)-1)*(nodes_y(2)-nodes_y(1)-1)), hwait);
    end;
    end;
    end;
    if ~isempty(hwait)
        close(hwait);
    end;
    disp(sprintf('File %s loaded', files(n).name));    
end;