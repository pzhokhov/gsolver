function refractive_index_raw_write(filename, lambda, n)
% wants lambda/n in order from high wavelengths to low wavelengths! In SI
% units.
fid = fopen(filename, 'wb');
fwrite(fid,3,'int'); 
N = length(lambda);
lambda = lambda(:);
if (~isreal(n))
    n_ = [real(n(:))'; imag(n(:))']; n_ = n_(:);
else
    n_ = [n(:)'; zeros(size(n(:)))']; n_ = n_(:);
end;
fwrite(fid, N, 'int');
const_SI; omega = 2*pi*SI.c./lambda;
fwrite(fid, omega, 'double');
Nw = fwrite(fid, n_, 'double'); 
if (Nw ~= 2*N) error('Error writing file!');end;
fclose(fid);

