function [W,L] = ion_aniso_(gamma, T, delta)

const_SI; 
lambda = 800e-9; 

Nt = 8192*2; 
t = -24*pi + (0:(Nt-1))*48*pi/Nt;  t=t(:);
t = t*1e-15/SI.atomic_time;
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


%w = 2*pi*SI.c/lambda/SI.atomic_frequency; 
w = 2e15/SI.atomic_frequency;


mx = 0.9;
my = 0.6; 

delta = delta*w; 


A0 = sqrt(mx*delta)./gamma;

E = real(exp(1i*w*t));  %electric field
fE = fft(E);  fE(abs(fE)<1e-3*max(abs(fE)))=0; 
E = ifft(fE);

A = cumsum(E)*dt; %vector potential
E = E.*A0./max(A); E0 = max(E); 
A = A./max(A).*A0; 


Ai = cumsum(A)*dt;
A2i=cumsum(A.^2)*dt; 

L = dt*cumsum(E./delta./(1+A.^2./delta)).*exp(-1i*delta*t);

W = abs(L).^2;