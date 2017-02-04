function W = calculate_plasma_keldysh(gamma, T, theta, delta)

const_SI; 
lambda = 800e-9; 

t = -30 + (0:8191)*60/8192; 
t = t*1e-15;
T = T*1e-15;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda; 


mx = 0.9;
my = 0.6;    

m = (cos(theta).^2./mx + sin(theta).^2./my)^(-1);


delta = delta*w*SI.h_;


E0 = sqrt(mx*SI.m*delta)*w/gamma/SI.e;
Ifactor = 2*sqrt(SI.mu0/SI.epsilon0);

E = E0*exp(-(t/T).^2);  %electric field amplitude
I = E.^2/Ifactor;

for i=1:length(I); 
    W(i) = keldysh_rate_ln(lambda, 1, delta, log(I(i)), m);
end;

W = cumsum(W)*(t(2)-t(1))*1e-15;


    