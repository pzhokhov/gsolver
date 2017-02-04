function n = refractive_index_raw(filename, omega)

fid = fopen(filename, 'rb');
type = fread(fid,1,'int'); 
if (type ~= 3) error('Incorrect type! Type is %d (must be 3)',type); end;
N = fread(fid, 1, 'int');
w = fread(fid, N, 'double');
n_ = fread(fid, 2*N, 'double'); n_=complex(n_(1:2:end),n_(2:2:end)); 

n = interp1(w, n_, omega);

fclose(fid);