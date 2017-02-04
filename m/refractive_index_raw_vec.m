function [no,ne] = refractive_index_raw_vec(filename, omega)

fid = fopen(filename, 'rb');
type = fread(fid,1,'int'); 
if (type ~= 6) error('Incorrect type! Type is %d (must be 6)',type); end;
N = fread(fid, 1, 'int');
w = fread(fid, N, 'double');
n_ = fread(fid, 4*N, 'double'); 
no_=complex(n_(1:4:end),n_(2:4:end)); 
ne_=complex(n_(3:4:end),n_(4:4:end));

ne = interp1(w, ne_, omega);
no = interp1(w, no_, omega);
