function [A,xnet, znet] = load_exdump(filename)

fid = fopen(filename, 'rb');
if (fid == -1) 
    error(sprintf('Unable to open file %s', filename));
end;

A = [];
znet = [];
N = fread(fid, 1, 'int');
xnet = fread(fid, N, 'double');
blockN = fread(fid, 1, 'int');

while (~feof(fid))
    z = fread(fid, 1,   'double');
    P = fread(fid, N, 'double');
    if (length(P) ~= N)
        break;
    end;
    A = [A P];
    znet = [znet z];
end;


fclose(fid);