function [d,t,E, c1, c2, c3] = atto_3hconv(cep, F, wup)

 Nt = 4096; 
 tmin = -200;
 tmax = 200;
 
t = tmin + (tmax-tmin)*(1:Nt)/Nt;
T = 3.5;
const_SI;
w0 = 2*pi*SI.c/800e-9/1e15;
windowlen = 2.5; 
%  
 %E = F*real(exp(-(t/T).^2+1i*t+1i*  cep));

 load('attotrans_data.mat'); 
 E = interp1(attotrans_data(:,1)*w0, attotrans_data(:,2), t, 'cubic', 0); 
 E = E./max(abs(E(:)))*F;  %
 E = E.*exp(-(t/max(attotrans_data(:,1))/w0/0.9).^8);
fE = fft(E); 
fE(1:(length(fE)/2))            = fE(1:(length(fE)/2))           *exp( 1i*cep); 
fE((length(fE)/2+1):length(fE)) = fE((length(fE)/2+1):length(fE))*exp(-1i*cep); 
fE(abs(fE)<1e-3*max(abs(fE)))=0;

E = ifft(fE);


Esim = E;
tsim = t; 

dt = 1; 
Q2 = 300; 

alpha = 0;
beta = 1;
w12 = 6.47;
tsim_ = tsim-min(tsim); 
R3 = zeros(size(tsim)); 
R3 = exp(-tsim*w12/Q2).*sin(w12*(tsim)).*(tsim>0);
%R((tsim-min(tsim)-t0) < 0) = 0;
%R3 = R3./sum(R3); 

R2 = ones(size(tsim_)).*(tsim>0);

R1 = zeros(size(tsim));
R1 = exp(-tsim*w12/Q2).*cos(w12*(tsim_)).*(tsim>0);

c3 = conv(Esim,     R3, 'same'); %c3=c3.*exp(-(t/max(attotrans_data(:,1))/w0/windowlen).^8);
c2 = conv(Esim.*c3, R2, 'same'); c2=c2.*exp(-(t/max(t)/0.8).^8);
c1 = conv(Esim.*c2, R1, 'same'); %c1=c1.*exp(-(t/max(attotrans_data(:,1))/w0/windowlen).^8);


d = alpha.*E.^3 + beta.*c1; 