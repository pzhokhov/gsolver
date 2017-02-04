function [E,t] = atto_E_peaks(h, I, T, W)

Nt = 8192;

tmin = -25;
tmax = 150;
t = tmin + (tmax-tmin)*(1:Nt)/Nt; t=reshape(t,length(t),1);
const_SI;

Emax = sqrt(2*I*sqrt(SI.mu0/SI.epsilon0))/SI.atomic_field;


Nt = length(t); 

if (nargin <= 2)
    T = 1.3488; %period (in fs)%
end;

if (nargin <= 3)
   W = 0.05*T;  %peak width
end;



t0 = 0; 
E = zeros(Nt, size(h,2)); 

dt = t(2)-t(1);

flux = Emax.^2*size(h,1)*(3.2439)^2*8;


for j=1:size(h,2); 
for i=1:size(h,1); 
    E(:,j) = E(:,j) + h(i,j)*exp(-((t-t0-(i-1)*T/2)./(W)).^2)*(2*mod(i,2)-1);
end;
E(:,j) = E(:,j)./sqrt(sum(abs(E(:,j))).^2*dt)*sqrt(flux);
%E(:,j)=E(:,j)*Emax./max(E(:,j));
fEj = (fft(E(:,j))); 
fEj(abs(fEj).^2<1e-3*max(abs(fEj)).^2)=0;
E(:,j)  = ifft(fEj); 
end;
