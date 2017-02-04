const_SI;
suffix='1e18_12f'; 
ftype = 'float'; 
I = 1e18; 

load Data_20121221_AllCEP_ReDone_20130408.mat
cep = Root.CEP(:).';
enu = Root.Energy(:);
eM  = squeeze(Root.Matrix(1,:,:)).';
Ncep = length(cep); 


[E, t]=atto_Ecep(cep*pi, I); 
Nt = size(E,1); 

w = cfreq(t); 
nu = SI.h_*w*1e15/SI.e; 
S = interp1(enu, eM, nu, 'linear', 0); 
S = fftshift(S,1); 
S = S./max(S(:)); 

fE = fftshift(fft(E),1); 
fEc = fE; fEc(w<0,:)=0;
Ec = ifft(ifftshift(fEc,1)); 

EcTHG = Ec.^3;
EcSPM = abs(Ec).^2.*Ec; 

fEcTHG = fft(EcTHG); 
fEcSPM = fft(EcSPM); 

fidE =fopen(sprintf('E_%s.bin', suffix), 'wb'); 
fwrite(fidE, Nt, 'int'); 
fwrite(fidE, Ncep, 'int');
fwrite(fidE, min(t), ftype); 
fwrite(fidE, max(t), ftype);
fwrite(fidE, E(:),   ftype); 
fclose(fidE); 

fidS =fopen(sprintf('S_%s.bin', suffix), 'wb'); 
fwrite(fidS, S(:),   ftype); 
fclose(fidS); 

fidTHG =fopen(sprintf('dTHG_%s.bin', suffix), 'wb'); 
T = [fEcTHG(1:2:end).'; fEcTHG(2:2:end).']; 
fwrite(fidTHG, T(:), ftype);
fclose(fidTHG); 


fidSPM =fopen(sprintf('dSPM_%s.bin', suffix), 'wb'); 
T = [fEcSPM(1:2:end).'; fEcSPM(2:2:end).']; 
fwrite(fidSPM, T(:), ftype);
fclose(fidSPM); 



