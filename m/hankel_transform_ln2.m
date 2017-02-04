    function [g,x] = hankel_transform_ln2(f, Nf)

f = f(:);
N = length(f);
phi = zeros(2*N,1);

alpha = log(2);
k0 = (2*exp(alpha)+exp(2*alpha))/(1+exp(alpha))^2/(1-exp(-2*alpha));
x0 = (1+exp(alpha))*exp(-alpha*N)/2;
x  = x0*exp(alpha*(0:(N-1))).';

phi_mult(1)    = k0*exp(alpha*(1-N));
phi_mult(2:(N))= exp(alpha.*(-N+(2:N)));

phi(1:(N)) = ([-diff(f(1:N)); f(N)].').*phi_mult;

j1  = (besselj(1,2*pi*Nf*x0*exp(alpha*((1-N):N)))).';
j1f = ifft(j1);

s = fft(fft(phi).*j1f);
g = 1./(x).*s(1:N); 