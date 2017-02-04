function [E,t] = atto_Ecep(cep, I)

Nt = 4096;

tmin = -25;
tmax = 120;
t = tmin + (tmax-tmin)*(1:Nt)/Nt; t=reshape(t,length(t),1);
const_SI;

Emax = sqrt(2*I*sqrt(SI.mu0/SI.epsilon0));


 load('attotrans_data1.mat'); 
 E0 = interp1(attotrans_data(:,1), attotrans_data(:,2), t, 'cubic', 0); 
 E0 = E0./max(abs(E0(:))).*Emax;  
 
 fE0 = (fft(E0)); 
 fE0(abs(fE0).^2<1e-3*max(abs(fE0)).^2)=0;
 E0  = ifft(fE0); 
 
 E = zeros(Nt, length(cep)); 
 
for j=1:length(cep);  E(:,j) = imag(hilbert(E0)*exp(1i*(pi*0.5+cep(j)))); end; 
 