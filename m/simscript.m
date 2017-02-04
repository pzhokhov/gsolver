function [d,t,E] = atto_3hconv(cep, F, wup)

 Nt = 4096; 
 tmin = -200;
 tmax = 200;
 
t = tmin + (tmax-tmin)*(1:Nt)/Nt;
T = 6;
const_SI;
w0 = 2*pi*SI.c/800e-9/1e15
%  
E = F*real(exp(-(t/T).^2+1i*t+1i*  cep));


Esim = E;
tsim = t; 

dt = 1; 
Q2 = 50; 

alpha = 0;
beta = 1;
w12 = 6.47
tsim_ = tsim-min(tsim); 
R3 = zeros(size(tsim)); 
R3 = exp(-tsim*w12/Q2).*sin(w12*(tsim)).*(tsim>0);
%R((tsim-min(tsim)-t0) < 0) = 0;
%R3 = R3./sum(R3); 

R2 = ones(size(tsim_)).*(tsim>0);

R1 = zeros(size(tsim));
R1 = exp(-tsim*w12/Q2).*cos(w12*(tsim_)).*(tsim>0)

c3 = conv(Esim,     R3, 'same'); 
c2 = conv(Esim.*c3, R2, 'same'); 
c1 = conv(Esim.*c2, R1, 'same'); 


d = alpha.*E.^3 + beta.*c3; 