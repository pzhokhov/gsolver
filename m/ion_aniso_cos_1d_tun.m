function [Wtun, t, E] = ion_aniso_cos(gamma, T, delta, cep, alpha)

if nargin <= 3
 cep = 0;
end;

if (nargin <= 4)
    alpha = 1; 
end; 

const_SI; 
lambda = 800e-9; 

Nt = 65536; 
t = -80 + (0:(Nt-1))*160/Nt;  t=t(:);
t = t*1e-15/SI.atomic_time;
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 1;
my = 1; 

delta = delta*w; 

alpha = 1; 

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


Itun = E.^2.*exp(-4*delta./abs(E)).*(pi*abs(E)).^(3/2)./sqrt(alpha); 
Wtun = cumsum(Itun)*dt;

return;

