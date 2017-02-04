function [W,t,E] = ion_1d_av(gamma, T, delta, cep, alpha)

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
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 1;
my = 1; 

delta = delta*w; 
 

A0 = 1./gamma;

E = real(exp(1i*w*t - (t/T).^2+1i*cep));  %electric field
%fE = fft(E);  fE(abs(fE)<1e-3*max(abs(fE)))=0; 
%E = ifft(fE);

A = cumsum(E)*dt; %vector potential

E0 = A0./w; 
E = E.*E0; E0 = max(E); 
A = A.*A0; 


% pmaxx = pi;
% px = -pmaxx + 2*pmaxx*(0:2047)./2048;
% 
% L = zeros(Nt, length(px));
% 
% Xi = t*ones(1,Nt); Xi=-Xi+Xi.';
% 
% E1E2 = E*E.';
% 
% 
% S = (2*pi).^3.*ion_Pfunction_3d(alpha.*delta.*t, A*cos(theta), A*sin(theta), zeros(size(A)));
% I = S.*exp(-1i.*delta.*(1+alpha).*Xi).*E1E2; 
% 
% W = zeros(size(t)); W(1)=I(1,1); 
% for nt=2:Nt; 
%     W(nt) = W(nt-1)+real(sum(I(nt,1:nt)))+real(sum(I(1:(nt-1),nt)));
% end; 
% W = W*dt.^2;
% return;


cA = cumsum(cos(A))*dt*alpha*delta;
sA = cumsum(sin(A))*dt*alpha*delta;
eA = cumsum(exp(1i*A))*dt*alpha*delta;
I = zeros(Nt,1);
for nt2 = 1:Nt; 
    C = cA(nt2) - cA;
    S = sA(nt2) - sA;
    
    Phi = sqrt(C.^2 + S.^2); 
    Phi = abs(eA(nt2)-eA);
    dPhi = Phi + alpha*delta*dt*((1:Nt).'-nt2); 
    
    %P = besselj(0, sqrt(C.^2 + S.^2));
    %P0 = besselj(0, alpha*delta*dt*((1:Nt).'-nt2));  
    %P1 = P - P0;
    
    %P =  P0;
    
    %K = cos(delta*(1+alpha)*(t(nt2)-t));
    %I(nt2) = 2*sum(P(1:nt2-1).*K(1:nt2-1).*E(1:nt2-1)) + P(nt2).*K(nt2).*E(nt2);
    
    P_ = real(exp(dPhi*(1i*pi/2)));
    %I(nt2) = 2*sum(P(1:nt2-1).*K(1:nt2-1).*E(1:nt2-1)) + P(nt2).*K(nt2).*E(nt2);
    I(nt2) = 2*sum(P_(1:nt2-1).*E(1:nt2-1)) + P_(nt2).*E(nt2);
end;

W = (2*pi).*cumsum(E.*I).*dt.^2;

return;

