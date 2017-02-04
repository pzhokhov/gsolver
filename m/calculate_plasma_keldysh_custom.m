function [W] = calculate_plasma_keldysh_a0(t_at, E, delta, alpha, d, lambda)

const_SI; 

if (nargin < 5)
 lambda = 800e-9; 
end;

w = 2*pi*SI.c/lambda;


Nt = length(t_at); 

t = t_at*(SI.h_/delta/SI.e);

m = SI.h_^2/d.^2/alpha/delta/SI.e/SI.m; 

E_SI = E*delta/d;

Ifactor = 2*sqrt(SI.mu0/SI.epsilon0);

I = (abs(E_SI)).^2/Ifactor;

delta_SI = delta.*SI.e;


parfor i=1:length(I); 
    W(i) = keldysh_rate_ln(lambda, 1, delta_SI, log(I(i)), m);
end;

W = cumtrapz(W)*(t(2)-t(1));


