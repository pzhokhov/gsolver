function [W,t,E, jc, jpa, jpa_] = ion_cos_1d_(w, A0, T, cep, alpha)

if nargin <= 3
 cep = 0;
end;

if nargin <= 4
  alpha = 1; 
end;

const_SI; 


Nt = 8192; 

T = T*2*pi; 
t = -6*T/w + 12*T/w*(0:(Nt-1))/Nt;
dt = t(2)-t(1); 
delta = 1; 


E = A0.*w.*real(exp(1i*w*t - (t*w/T).^2+1i*cep));  %electric field
E = E(:); 
t = t(:);
A = cumtrapz(E)*dt; %vector potential

cA = cumtrapz(cos(A))*dt*alpha*delta;
sA = cumtrapz(sin(A))*dt*alpha*delta;


jpa = zeros(Nt,1); 
Ijc_c = zeros(Nt,1); 
Ijc_s = zeros(Nt,1); 

I = zeros(Nt,1);
parfor nt2 = 1:Nt; 
    C = cA(nt2) - cA;
    S = sA(nt2) - sA;
    
    P0 = besselj(0, sqrt(C.^2 + S.^2));
    P1 = besselj(1, sqrt(C.^2 + S.^2)); 
    argphi = atan2(S,C); 
    
    
    K = exp(1i*delta*(1+alpha)*(t(nt2)-t));
    
    I(nt2) = 2*trapz(P0(1:nt2).*real(K(1:nt2)).*E(1:nt2));
    
    Ijc_c(nt2) = 2*trapz(P1(1:nt2).*imag(K(1:nt2)).*E(1:nt2).*cos(argphi(1:nt2)));
    Ijc_s(nt2) = 2*trapz(P1(1:nt2).*imag(K(1:nt2)).*E(1:nt2).*sin(argphi(1:nt2)));
    
    jpa(nt2) = 4*pi*dt*trapz(E(1:nt2).*delta.*(real(K(1:nt2)).*(1+alpha).*P0(1:nt2) - alpha.*P1(1:nt2).*cos(argphi(1:nt2)).*imag(K(1:nt2))));
    %jpa(nt2) = 2*dt*(trapz(E(1:nt2).*real(K(1:nt2)).*(delta.*(1+alpha).*P0(1:nt2))));
    
 end;

W  = (2*pi).*cumtrapz(E.*I).*dt.^2;
jc = (2*pi)./2.*delta.*alpha.*dt.^2*(cumtrapz(E.*Ijc_s(:)).*cos(A) - cumtrapz(E.*Ijc_c(:)).*sin(A));

jpa_ = 2*pi*I.*dt.*delta.*(1+alpha);
return;

