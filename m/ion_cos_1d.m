function [W,t,E, sigmac, sigmapa, jc, jpa] = ion_cos_1d(a0, T, w, cep, alpha)

if nargin <= 3
 cep = 0;
end;

if nargin <= 4
  alpha = 1; 
end;


Nt = 2048; 

delta = 1; %delta*w; 

T=T/sqrt(log(4))*(2*pi/w); 

t = -6*T + (0:(Nt-1))/Nt*12*T;
t = t.';

dt = t(2)-t(1);

A0 = a0; 

E = A0.*w.*real(exp(1i*w*t - (t/T).^2+1i*cep));  %electric field

A = cumsum(E)*dt; %vector potential

cA = cumsum(cos(A))*dt*alpha*delta;
sA = cumsum(sin(A))*dt*alpha*delta;


jpa = zeros(Nt,1); 
Ijc_c = zeros(Nt,1); 
Ijc_s = zeros(Nt,1); 

I = zeros(Nt,1);

sigmac = zeros(Nt, Nt); 
sigmapa = zeros(Nt, Nt); 

tic_gp; 
for nt2 = 1:Nt; 
    C = cA(nt2) - cA;
    S = sA(nt2) - sA;
    
    P0 = besselj(0, sqrt(C.^2 + S.^2));
    P1 = besselj(1, sqrt(C.^2 + S.^2)); 
    argphi = atan2(S,C); 
    
    
    K = exp(1i*delta*(1+alpha)*(t(nt2)-t));
    
    I(nt2) = 2*trapz(P0(1:nt2).*real(K(1:nt2)).*E(1:nt2));
    
    Ijc_c(nt2) = 2*trapz(P1(1:nt2).*imag(K(1:nt2)).*E(1:nt2).*cos(argphi(1:nt2)));
    Ijc_s(nt2) = 2*trapz(P1(1:nt2).*imag(K(1:nt2)).*E(1:nt2).*sin(argphi(1:nt2)));
    
    
    
    jpa(nt2) = 4*pi*dt*trapz(E(1:nt2).*real(K(1:nt2)).*(delta.*(1+alpha).*P0(1:nt2) - delta.*alpha.*P1(1:nt2).*imag(K(1:nt2)).*cos(argphi(1:nt2))));
    
    sigmapa(nt2, 1:nt2) = 4*pi.*real(K(1:nt2)).*(delta.*(1+alpha).*P0(1:nt2) - delta.*alpha.*P1(1:nt2).*imag(K(1:nt2)).*cos(argphi(1:nt2)));
    %jpa(nt2) = 2*dt*(trapz(E(1:nt2).*real(K(1:nt2)).*(delta.*(1+alpha).*P0(1:nt2))));
        
    toc_gp(nt2/Nt);
 end;
sigmac = (2*pi)./2.*delta.*dt.*(cos(A)*Ijc_s.' - sin(A)*Ijc_c.');
for nt=1:Nt;
    sigmac(nt,nt+1:end) = 0;
end;
W  = (2*pi).*cumtrapz(E.*I).*dt.^2;
jc = (2*pi)./2.*delta.*dt.^2*(cumtrapz(E.*Ijc_s(:)).*cos(A) - cumtrapz(E.*Ijc_c(:)).*sin(A));



%jpa = I;
t = w*t/2/pi; 
return;

