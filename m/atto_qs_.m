
function [d, t, E, rhop, rhopp] = atto_qs(cep, I, e12, d12, Q2, al, Efun)

Nt = 8192;

tmin = -25;
tmax = 300;
t = tmin + (tmax-tmin)*(1:Nt)/Nt; t=reshape(t,length(t),1);
const_SI;

%w0 = 2*pi*SI.c/800e-9/1e15;
w12 = e12.*SI.e/SI.atomic_energy;


%d12 = (3*SI.h_.*f12/2/SI.m./w12/1e15).^(0.5)*SI.e;
Emax = sqrt(2*I*sqrt(SI.mu0/SI.epsilon0))/SI.atomic_field;

[E,t] = Efun(cep, I); 
 
Emax = sqrt(2*I*sqrt(SI.mu0/SI.epsilon0))/SI.atomic_field;
E = E./max(abs(E)).*Emax;
 

nlev = numel(w12); 
Nt = length(t); 
w12 = reshape(w12, 1, nlev); 
d12 = reshape(d12, 1, nlev); 

Omega = d12; %%*Emax;

gamma12=w12./Q2; 

if (nargin < 6)
    Efun = atto_Ecep; 
end;

%  load('attotrans_data1.mat'); 
%  E0 = interp1(attotrans_data(:,1), attotrans_data(:,2), t, 'cubic', 0); 
%  E0 = E0./max(abs(E0(:)));  
%  
%  fE0 = (fft(E0)); 
%  fE0(abs(fE0).^2<1e-3*max(abs(fE0)).^2)=0;
%  E0  = ifft(fE0); 
%  
%  E = zeros(Nt, length(cep)); 

 

 
 t = t*1e-15/SI.atomic_time; %convert t into the atomic units of time
 
% E = (ifft(fE));
d = zeros(size(E));
dt = t(2)-t(1); 
di = zeros(Nt, nlev, length(cep));
dn = zeros(Nt, length(cep));
for j=1:length(cep); 
 cE = E(:,j);
 %cE = imag(hilbert(E0)*exp(1i*(pi*0.5+cep(j))));
 %E(:,j)=cE; 
 %Edot = [0; cE(3:end)-cE(1:end-2); 0]/2;
 Edot = gradient(cE); 
 
 weff = sqrt   ((ones(Nt,1)*w12).^2 + 4*(cE*(Omega)).^2); 
 weff = reshape( weff, length(weff), nlev);  
 
 rhopp = zeros(Nt, nlev); %rhop=rhopp; 

 for i=1:Nt;
   theta_t = (t<=t(i)); 
   Theta_t = theta_t*ones(1,nlev);
   
   Phi = (ones(Nt,1)*trapz(weff.*Theta_t,1) - cumtrapz(weff.*(Theta_t))).*dt;
   R = (Theta_t.*sin(Phi)).*exp(-(t(i)-t)*gamma12);
   rhopp(i,:) = (trapz(R.*(Edot*Omega),1));
 end; 
 
rhopp = rhopp./weff; 
rhop = w12.*cumtrapz(rhopp)*dt; 
d = 2*rhop*d12;

%dij = cumtrapz(rhopp./weff,1).*dt.*(exp(-(t/max(t)/0.8).^10)*(al.*d12));
%dij = cumtrapz(rhopp./weff,1).*dt.*(exp(-(t/max(t)/0.8).^10)*(al.*d12));
%dn(:,j) = cumtrapz(rhopp./w12.*cE*Omega)*dt;
%dij = cumsum(rhopp./weff,1).*dt.*(ones(Nt,1)*(al.*d12)); 
%d(:,j)=squeeze(sum(dij, 2));
%di(:,:,j)=dij;
end; 
%di = permute(di, [1,3,2]); 