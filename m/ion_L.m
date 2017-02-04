function [W,L,t] = ion_L(gamma, T, delta,  cep, alph);


if nargin <= 3
 cep = 0;
end;

if (nargin <= 4)
 alph = 1; 
end;

const_SI; 
lambda = 800e-9; 

Nt = 2^18; 
t = -40 + (0:(Nt-1))*80/Nt;  t=t(:);
t = t*1e-15/SI.atomic_time;
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 1;
my = 1; 

delta = delta*w; 


A0 = 1./gamma./sqrt(alph);

E = real(exp(1i*w*t - (t/T).^2+1i*cep));  %electric field
%fE = fft(E);  fE(abs(fE)<1e-3*max(abs(fE)))=0; 
%E = ifft(fE);

A = cumsum(E)*dt; %vector potential

E0 = A0./w; 
E = E.*E0; E0 = max(E); 
A = A.*A0; 



%pmaxx = sqrt(mx*delta); 
Np = 2048;
px = -pi + (0:(Np-1))*2*pi/Np;

L = zeros(Nt, length(px));
for nx = 1:Np;
    %en = delta + (px(nx)+A).^2;
    en = delta*(1+ alph*(px(nx)+A).^2);
   
    K = exp(-1i*cumsum(en)*dt);
    L(:, nx) = dt*cumsum(K(:).*E(:)); 
end;
      

W = sum(abs(L).^2, 2)*dt; 
