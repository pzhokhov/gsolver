function [W,W2, L,E, t] = ion_aniso_cos(gamma, T, delta, theta)

if nargin <= 3
 theta = 0;
end;

const_SI; 
lambda = 800e-9; 

Nt = 2048; 
t = -40 + (0:(Nt-1))*80/Nt;  t=t(:);
t = t*1e-15/SI.atomic_time;
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 0.9;
my = 0.6; 

delta = delta*w; 

alpha = 0.2; 

A0 = 1./gamma./sqrt(alpha);

E = real(exp(1i*w*t - (t/T).^2+1i*pi/2));  %electric field
fE = fft(E);  fE(abs(fE)<1e-3*max(abs(fE)))=0; 
E = ifft(fE);

A = cumsum(E)*dt; %vector potential

E = E.*A0./max(A); E0 = max(E); 
A = A./max(A).*A0; 


pmaxx = pi;
px = -pmaxx + 2*pmaxx*(0:2047)./2048;

L = zeros(Nt, length(px));

Xi = t*ones(1,Nt); Xi=-Xi+Xi.';

E1E2 = E*E.';




% Ic_ = cumsum(cos(A)); 
% Is_ = cumsum(sin(A)); 
% 
% Ic = Ic_(:)*ones(1,Nt); Ic = -Ic + Ic.';
% Is = Is_(:)*ones(1,Nt); Is = -Is + Is.';
% Ix = alpha.*delta.*sqrt(Ic.^2 + Is.^2);

%S = besselsum(alpha.*delta.*Ic, alpha.*delta.*Is); 
S = (2*pi).^3.*ion_Pfunction_3d(alpha.*delta.*t, A*cos(theta), A*sin(theta), zeros(size(A)));
I = S.*exp(-1i.*delta.*(1+alpha).*Xi).*E1E2; 

W = zeros(size(t)); W(1)=I(1,1); 
for nt=2:Nt; 
    W(nt) = W(nt-1)+real(sum(I(nt,1:nt)))+real(sum(I(1:(nt-1),nt)));
end; 
W = W*dt.^2;
return;

parfor nx = 1:length(px);
    en = delta.*(1+alpha - alpha.*cos(px(nx)+A));
    %ien = delta.*t_ + A2i + (px(nx)).^2.*t_ ;
    K = exp(-1i*cumsum(en).*dt)
    L(:, nx) = dt*cumsum(K(:).*E(:)); 
end;
W2 = sum(abs(L).^2, 2)*(px(2)-px(1)); 

return; 

