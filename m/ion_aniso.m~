    
function [W, L, t, A, E, px, py] = ion_aniso(gamma, T, theta, delta)

const_SI; 
lambda = 800e-9; 

t = -30 + (0:8191)*60/8192; 
t = t*1e-15/SI.atomic_time;
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 0.9;
my = 0.6; 

delta = delta*w; 


A0 = sqrt(mx*delta)./gamma;

%T = 10;
%E0 = 10;



E = real(exp(1i*w*t - (t/T).^2));  %electric field
fE = fft(E);  fE(abs(fE)<1e-3*max(abs(fE)))=0; 
E = ifft(fE);

A = cumsum(E)*dt; %vector potential
E = E.*A0./max(A); E0 = max(E); 
A = A./max(A).*A0; 



if (nargin <= 2)
  ex_ = 1; ey_ = 1; %unnormalized polarization vector
else
  ex_ = cos(theta); ey_ = sin(theta); 
end;

ex =  ex_./sqrt(ex_.^2 + ey_.^2); 
ey =  ey_./sqrt(ex_.^2 + ey_.^2); 


L = zeros(length(t), length(px)); 
clear en; 

for nx = 1 : length(px)
    
    for nt = 1 : length(t); 
        en(nt) = delta + (px(nx)+A(nt)).^2;
    end; 
    K = exp(-1i*cumsum(en)*dt)./en; 
    
    L(:, nx) = dt*cumsum(K(:).*E(:)./E0);
end;

dpx = px(2)-px(1); 

W = sum(abs(L(:,:)).^2,2)).*dpx;