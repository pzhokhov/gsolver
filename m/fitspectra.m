    
function [al, mS, Ximin] = fitspectra(nu, S0, S); 

xstart = [ones(1, size(S,2))];
[x0, Ximin] = fminsearch(@(x)(Xi2(x, S0, S)), xstart, optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4)); 

Nx = size(x0, 2); 
a0 = x0(1:Nx); 
phi0 = zeros(size(x0)); %x0(Nx/2+1:Nx);
al = a0.*exp(1i*phi0); 

mS = abs(sum(S.*(ones(size(S,1),1)*al),2)).^2;
Ximin

function f=Xi2(x, S0, S); 

Nx = size(x,2); 
Nw = size(S,1); 
a = x; 
phi = zeros(size(a)); %x(Nx/2+1:Nx);
mS = abs(sum(S.*(ones(Nw,1)*(a.*exp(1i*phi))),2)).^2;
f = sum(abs(S0 - mS).^2);
return; 




