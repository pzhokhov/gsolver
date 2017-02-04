function H = hankel_transform_quad(F);

F = F(:);

N = size(F,1);

for l=1:N;
  S = besselj(0, 2*pi*(l-1)*(0:(N-1))/N).'.*(0:(N-1)).'.*F; 
  H(l) = trapz(S).*2*pi/N/1.3092/1.0871;
end;

