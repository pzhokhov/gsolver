function R = phaseshift(X, phi, dim)

if (nargin < 3)
    dim = 1;
end;
if (dim ~= 1)
     Ndim = length(size(X)); 
     D=[dim, 2:(dim-1), 1, dim+1:Ndim];
     X = permute(X, D); 
end;

N = size(X,1);

fX = fft(X); 
fX_ = zeros(size(X)); 

fX_(1:N/2,:)=fX(1:N/2,:)*exp(-1i*phi); 

fX_(N/2+1:end,:) = fX(N/2+1:end,:)*exp(1i*phi); 

R = ifft(fX_); 

if (dim ~= 1)
    R = permute(R, D);
end;
