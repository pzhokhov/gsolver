function [W,S] = ion_cos_av(gamma, T, delta, theta, cep)

if nargin <= 3
 theta = 0;
end;

if nargin <= 4
 cep = 0;
end;

const_SI; 
lambda = 800e-9; 

Nt = 8192; 
t = -80 + (0:(Nt-1))*160/Nt;  t=t(:);
t = t*1e-15/SI.atomic_time;
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 1;
my = 1; 

delta = delta*w; 

alpha = 1; 

A0 = 1./gamma./sqrt(alpha);

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


Ax = A*cos(theta); 
Ay = A*sin(theta);
Az = zeros(size(A)); 

cAppp = cumsum(cos(Ax+Ay+Az))*dt/4*alpha*delta;
sAppp = cumsum(sin(Ax+Ay+Az))*dt/4*alpha*delta;

cAppm = cumsum(cos(Ax+Ay-Az))*dt/4*alpha*delta;
sAppm = cumsum(sin(Ax+Ay-Az))*dt/4*alpha*delta;

cApmp = cumsum(cos(Ax-Ay+Az))*dt/4*alpha*delta;
sApmp = cumsum(sin(Ax-Ay+Az))*dt/4*alpha*delta;

cApmm = cumsum(cos(Ax-Ay-Az))*dt/4*alpha*delta;
sApmm = cumsum(sin(Ax-Ay-Az))*dt/4*alpha*delta;

I = zeros(Nt,1);
parfor nt2 = 1:Nt; 
    Cppp = cAppp(nt2) - cAppp;
    Sppp = sAppp(nt2) - sAppp;
    
    Cpmp = cApmp(nt2) - cApmp;
    Spmp = sApmp(nt2) - sApmp;
    
    Cppm = cAppm(nt2) - cAppm;
    Sppm = sAppm(nt2) - sAppm;
    
    Cpmm = cApmm(nt2) - cApmm;
    Spmm = sApmm(nt2) - sApmm;
    
    P = besselsum4(Cppp, Sppp, Cpmp, Spmp, Cppm, Sppm, Cpmm, Spmm); 
    K = cos(delta*(1+alpha)*(t(nt2)-t));
    I(nt2) = 2*sum(P(1:nt2-1).*K(1:nt2-1).*E(1:nt2-1)) + P(nt2).*K(nt2).*E(nt2);
end;
W = (2*pi).^3.*cumsum(E.*I).*dt.^2;

return;

parfor nx = 1:length(px);
    en = delta.*(1+alpha - alpha.*cos(px(nx)+A));
    %ien = delta.*t_ + A2i + (px(nx)).^2.*t_ ;
    K = exp(-1i*cumsum(en).*dt)
    L(:, nx) = dt*cumsum(K(:).*E(:)); 
end;
W2 = sum(abs(L).^2, 2)*(px(2)-px(1)); 

return; 

