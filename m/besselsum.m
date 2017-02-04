function S=besselsum(A,B,k);

S = zeros(size(A));
if (nargin<3)
  parfor m=0:9;
    M = m.*ones(size(A));
    T = besselj(M,A).*besselj(M,B).*((1i)^m); 
    if (m>0) T = T.*2; end;
    S = S+T;
  end;
  return; 
else
  parfor m=-9:9;
    M = m.*ones(size(A));
    T = besselj(M,A).*besselj(M+k,B).*((1i)^m); 
    S = S+T;
  end;
  return;
end;