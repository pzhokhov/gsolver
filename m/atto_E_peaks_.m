function [E,t, E0, t_at] = atto_E_peaks_(cep, I)

const_SI;

[E0,t] = atto_E_peaks([1;1], I); 
Nt = length(t);

 fE0 = (fft(E0)); 
 fE0(abs(fE0).^2<1e-3*max(abs(fE0)).^2)=0;
 E0  = ifft(fE0); 
 
 E = zeros(Nt, length(cep)); 
 
 t_at = t*1e-15/SI.atomic_time; %convert t into the atomic units of time
 
% E = (ifft(fE));

for j=1:length(cep); 
 cE = imag(hilbert(E0)*exp(1i*(pi*0.5+cep(j))));
 E(:,j)=cE; 
end; 
