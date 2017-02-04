function [d] = atto_h(cep, I, Efun)



    
const_SI;

 peaks = cep; 
 Ncep = size(peaks,2); 
 
 if (nargin <= 2)
     Efun = @(a,b)atto_E_peaks(a,b,8); 
 end;
 
 [E,t] = Efun(peaks, I);

 w = cfreq(t); 
 nu = w.*SI.h_./SI.e.*1e15; 
 
fE = fftshift(fft(E),1); 
fEp = fE; fEp(nu<0,:)=0;
Ep = ifft(ifftshift(fEp,1)); 

Ep3 = Ep.^3; 
Epspm = Ep.^2.*conj(Ep); 
Ep5 = Ep.^5; 

d = Epspm/5 + Ep3/10 + Ep5.*exp(-1i*pi/2)*5; 


