function [d, fd] = atto_qsr(cep, F, wup)

 Nt = 2048;
 tmin = -30;
 tmax = 50;
 
t = tmin + (tmax-tmin)*(1:Nt)/Nt; t=reshape(t,length(t),1);
T = 4*pi;
const_SI;
w0 = 2*pi*SI.c/800e-9/1e15;
windowlen = 2.5; 

%  


P = (abs(t/T)<1/2); 

S=0;
Q = 200;
z=F.^2/2/wup; 
weff = wup + 2*z;
Gamma = wup/Q; 

for k=-3:3; 
for n=-3:3; 
        S = S + (besselj(n, z)-besselj(n-1, z)).*besselj(k-n, z)./(-1i*(weff + (2*n-1))+Gamma).*exp(1i.*(2*k-1).*(t+cep)); 
        S = S + besselj(k,z)*((besselj(n-1,z)+besselj(n,z))/1i/T + (besselj(n-1,z)-besselj(n,z))/(-1i*(weff + (2*n-1))+Gamma)).*(exp(1i.*weff.*(t+T/2)+2i*k.*(t+cep) - i*(2*n-1).*(-T/2+cep)));
        
end;
end;

d = P.*F/4.*S;
d = d + conj(d); 
fd = fftshift(fft(d),1);