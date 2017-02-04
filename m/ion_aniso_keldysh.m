function [W] = ion_aniso_keldysh(gamma, T, delta, theta)

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


mx = 1;
my = 0.6; 

delta = delta*w; 


A0 = sqrt(mx*delta)./gamma;

Ec = A0.*w.*(exp(1i*w*t - (t/T).^2)); 
alpha = 0.2; 


dW = zeros(Nt, 1); 
gammat = gamma./abs(Ec(:)).*max(abs(Ec));
for nt=1:Nt;
    dW(nt) = keldysh_rate_au(lambda, delta/w, gammat(nt), 1); 
end;

W = cumsum(dW)*dt;
