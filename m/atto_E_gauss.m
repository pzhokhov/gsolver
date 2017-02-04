function [E,t, E0, t_at] = atto_E_gauss(cep, I)

const_SI;


Nt = 8192;

tmin = -25;
tmax = 150;
t = tmin + (tmax-tmin)*(1:Nt)/Nt; t=reshape(t,length(t),1);
const_SI;

Emax = sqrt(2*I*sqrt(SI.mu0/SI.epsilon0))/SI.atomic_field;


Nt = length(t); 

l0 = 420e-9; 
w0 = 2*pi*SI.c/l0/1e15; 
T = 0.7;

E0 = Emax.*real(exp(-(t/T).^2+1i*w0*t));


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
