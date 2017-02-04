    
function [d, al, dfit] = fitspectra_(ecep, enu, eM, t, xstart, I, nulev, Qs); 


dstart = xstart(1:((length(xstart)-2)/2));

[dfit, dfit_, E, t] = atto_qs(ecep*pi, 1e18, nulev, dstart, Qs, 1);
fdfit  = fftshift(fft(dfit),1); 
fdfit_ = fftshift(fft(dfit_),1); 
fE = fft(E); ; 
Nt = size(E,1); 
fE(Nt/2:Nt,:)=0; 
Ec = ifft(fE); 
fEc_THG = fftshift(fft(Ec.^3),1); 
fEc_SPM = fftshift(fft(abs(Ec).^2.*Ec),1); 

const_SI; 
w = cfreq(t); 
nu = w*SI.h_./SI.e*1e15; 

s = 6 < nu  & nu < 14; 
nstart = 13; 

M0  = interp1(enu, eM(:,:), nu, 'linear', 0);

M0 = M0./max(max(abs(M0))); 
S0 = M0(:,nstart);

%[alstart, mS, Ximin] = fitspectra(nu(s), S0(s), [squeeze(fdfit_(s,nstart,:))]);

%xstart = [dstart   alstart 0 0]; 

tic;
[x, ximin] = fminsearch(@(x)xi2(x, s, M0, ecep, I, nulev, Qs, fEc_SPM, fEc_THG), xstart, optimset('MaxFunEvals', 500, 'MaxIter', 500));


[dfit, dfit_, E, t] = atto_qs(ecep*pi, 1e18, nulev, dstart, Qs, 1);





function f=xi2(x, s, M, ecep, I, nulev, Qs, fEc_SPM, fEc_THG); 

Nx = size(x,2); 
f12 = x(1:(Nx/2-1)); al = x((Nx/2):(Nx-2)); 
B = x(end-1); 
C = x(end); 

%[dfit, ~, E] = atto_qs(ecep*pi, I, nulev, f12, Qs, al);

Nt = size(E,1); 
fdfit = fftshift(fft(dfit),1); fdfit(~s,:)=0; fdfit = fdfit./max(abs(fdfit(:)));

f = sum(sum((abs(M)-abs(fdfit+B*fEc_SPM+C*fEc_THG).^2).^2));
f 
x
toc;
return; 




