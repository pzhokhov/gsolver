function [A,xnet, znet] = load_exdumpold(filename)

fid = fopen(filename, 'rb');
if (fid == -1) 
    error(sprintf('Unable to open file %s', filename));
end;

A = [];
znet = [];
N = fread(fid, 1, 'int');
xnet = fread(fid, N, 'double');


while (~feof(fid))
    z = fread(fid, 1,   'double');
    P = fread(fid, 2*N, 'double');
    if (length(P) ~= 2*N)
        break;
    end;
    A = [A P(1:end/2)];
    znet = [znet z];
end;


fclose(fid);