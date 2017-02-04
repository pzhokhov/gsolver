function [d,E,t, R] = atto_qs__(cep, F, wup, dphi, Q)

 Nt = 4096;
 tmin = -50;
 tmax = 150;
 
t = tmin + (tmax-tmin)*(1:Nt)/Nt; t=reshape(t,length(t),1);
T = 4*pi;
const_SI;
w0 = 2*pi*SI.c/800e-9/1e15;
windowlen = 2.5; 

w = cfreq(t);
%  
% E = F*real(exp(-(t/T).^2+1i*t+1i*  cep));

%E = F.*cos(t + cep).*(abs(t/T)<1/2); 

 load('attotrans_data.mat'); 
     E = interp1(attotrans_data(:,1)*w0, attotrans_data(:,2), t, 'cubic', 0); 
     E = E./max(abs(E(:)))*F;  %
     E = E.*exp(-(t/max(attotrans_data(:,1))/w0/0.9).^8);
    fE = (fft(E)); 
    fE(1:(length(fE)/2))            = fE(1:(length(fE)/2))           *exp( 1i*cep); 
    fE((length(fE)/2+1):length(fE)) = 0; %fE((length(fE)/2+1):length(fE))*exp(-1i*cep); 
    fE(abs(fE)<1e-3*max(abs(fE)))=0;
    %fE(w<0)=0;

    E = (ifft(fE));


Esim = E;
tsim = t; 

dt = 1; 
Q2 = Q; 

alpha = 0;
beta = 1;
w12 = wup;
gamma12=w12/Q2; 
% tsim_ = tsim-min(tsim); 
% R3 = zeros(size(tsim)); 
% R3 = exp(-tsim*w12/Q2).*sin(w12*(tsim)).*(tsim>0);
% %R((tsim-min(tsim)-t0) < 0) = 0;
% %R3 = R3./sum(R3); 
% 
% R2 = ones(size(tsim_)).*(tsim>0);
% 
% R1 = zeros(size(tsim));
% R1 = exp(-tsim*w12/Q2).*cos(w12*(tsim_)).*(tsim>0);
%   
% c3 = conv(Esim,     R3, 'same'); %c3=c3.*exp(-(t/max(attotrans_data(:,1))/w0/windowlen).^8);
% c2 = conv(Esim.*c3, R2, 'same'); c2=c2.*exp(-(t/max(t)/0.8).^8);
% c1 = conv(Esim.*c2, R1, 'same'); %c1=c1.*exp(-(t/max(attotrans_data(:,1))/w0/windowlen).^8);


%weff = sqrt(w12.^2 + 4*Esim.^2); 
%weff = w12 + 2*Esim.^2./w12; weff=reshape(weff, length(weff),1); 
%Esimdot = [0; diff(Esim)];
%Esimdot = Esim;

Hs = abs(E).^2.*E + E.^3/3;
Hc = E.^5; 
d=zeros(length(t),1); 
Rs=zeros(length(t)); 
Rc=Rs; 
for i=1:length(t);
   phi  = (t(i)-t)*wup;
   phiu = phi + ones(size(t))*dphi;
   Rc(:,i) = (t(i)>t).*sum(exp(1i*phi  - (t(i)-t)*gamma12),2);
   Rs(:,i) = (t(i)>t).*sum(exp(1i*phiu - (t(i)-t)*gamma12),2);
   d(i) = -sum(Rc(:,i).*(Hc) + Rs(:,i).*Hs);
end;
