function [E, t, E_0, t_at] = atto_E_exp(cep, I)

Nt = 8192;

tmin = -25;
tmax = 25;
t = tmin + (tmax-tmin)*(1:Nt)/Nt; t=reshape(t,length(t),1);
const_SI;

Emax = sqrt(2*I*sqrt(SI.mu0/SI.epsilon0))/SI.atomic_field;


Nt = length(t); 

 load('attotrans_data1.mat'); 
 E0 = interp1(attotrans_data(:,1), attotrans_data(:,2), t, 'cubic', 0); 
 E0 = E0./max(abs(E0(:)));  
 
 fE0 = (fft(E0)); 
 fE0(abs(fE0).^2<1e-3*max(abs(fE0)).^2)=0;
 E0  = Emax.*ifft(fE0); 
 
 E = zeros(Nt, length(cep)); 
 
 t_at = t*1e-15/SI.atomic_time; %convert t into the atomic units of time
 
% E = (ifft(fE));

for j=1:length(cep); 
 cE = imag(hilbert(E0)*exp(1i*(pi*0.5+cep(j))));
 E(:,j)=cE; 
end; 