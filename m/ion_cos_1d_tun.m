function [W,S] = ion_cos_1d_tun(gamma, T, delta, cep, alpha)

if nargin <= 3
 cep = 0;
end;

if nargin <= 4
  alpha = 1; 
end;

const_SI; 
lambda = 800e-9; 

Nt = 16384; 
t = -80 + (0:(Nt-1))*160/Nt;  t=t(:);
t = t*1e-15/SI.atomic_time;
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 1;
my = 1; 

delta = delta*w; 
 

A0 = 1./gamma;

E = real(exp(1i*w*t - (t/T).^2+1i*cep));  %electric field
Ec = exp(1i*w*t - (t/T).^2+1i*cep); 
%fE = fft(E);  fE(abs(fE)<1e-3*max(abs(fE)))=0; 
%E = ifft(fE);

A = cumsum(E)*dt; %vector potential
Ac = cumsum(Ec).*dt.*A0; 

E0 = A0./w; 
E = E.*E0; E0 = max(E); 
A = A.*A0; 


deltat = delta.*(1 + alpha.*abs(Ac).^2./4); 
ionrate = 2*sqrt(2*pi)*abs(Ec).^2.*real(exp(1i*w.*t - 1i*deltat.*t + 1i.*alpha.*delta.*abs(Ac).^2./8./w.*sin(2.*w*t) - 1i*pi*4 -...
          2*deltat./w./sqrt(alpha)./abs(Ac))./2^(1/4)./alpha^(3/4)./delta./sqrt(w.*t - 1i.*sqrt(2)./sqrt(alpha)./abs(Ac)));  


W = cumsum(ionrate).*dt;

return;

