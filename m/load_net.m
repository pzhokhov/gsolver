function [net] = load_net(filename)

fid = fopen(filename, 'rb');
N = fread(fid, 1, 'int');
net = fread(fid, N, 'double');

fclose(fid);