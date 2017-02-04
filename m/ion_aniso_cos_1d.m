function [W,t,E, jc, jpa] = ion_aniso_cos_1d(gamma, T, delta, cep, alpha)

if nargin <= 3
 cep = 0;
end;

if nargin <= 4
  alpha = 1; 
end;

const_SI; 
lambda = 800e-9; 

Nt = 16384; 
t = -40 + (0:(Nt-1))*80/Nt;  t=t(:);
t = t*1e-15/SI.atomic_time;
%T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 1;
my = 1; 

delta = delta*w; 
 

A0 = 1./gamma/sqrt(alpha);

E = A0.*w.*real(exp(1i*w*t - (t*w/T).^2+1i*cep));  %electric field

A = cumsum(E)*dt; %vector potential

cA = cumsum(cos(A))*dt*alpha*delta;
sA = cumsum(sin(A))*dt*alpha*delta;

I = zeros(Nt,1);
parfor nt2 = 1:Nt; 
    C = cA(nt2) - cA;
    S = sA(nt2) - sA;
    
    P0 = besselj(0, sqrt(C.^2 + S.^2));
    P1 = besselj(1, sqrt(C.^2 + S.^2)); 
    argphi = atan2(S,C); 
    
    
    K = exp(1i*delta*(1+alpha)*(t(nt2)-t));
    
    I(nt2) = 2*sum(P0(1:nt2-1).*real(K(1:nt2-1)).*E(1:nt2-1)) + P0(nt2).*real(K(nt2)).*E(nt2);
    
    
end;

W  = (2*pi).*cumsum(E.*I).*dt.^2;
W1 = (2) 


return;

