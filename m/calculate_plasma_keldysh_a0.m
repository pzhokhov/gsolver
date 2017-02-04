function [W,t, cE, jck] = calculate_plasma_keldysh_a0(a0, T, delta, alpha, d, lambda, sigma)

const_SI; 

if (nargin < 5)
 lambda = 800e-9; 
end;

w = 2*pi*SI.c/lambda;

T = T*2*pi/sqrt(log(4)); 

Nt = 8192; 
t = -6*T/w + 12*T/w*(0:(Nt-1))/Nt;

m = SI.h_^2/d.^2/alpha/delta/SI.e/SI.m; 

Emax = (w*SI.h_.*a0./SI.e/d);
Emax_at = a0*w*SI.h_/delta/SI.e; 

cE =       Emax.* exp(-(w.*t/T).^2+ 1i*w.*t);
cE_at = Emax_at.* exp(-(w.*t/T).^2+ 1i*w.*t);


Ifactor = 2*sqrt(SI.mu0/SI.epsilon0);

I = (abs(cE)).^2/Ifactor;
deltaSI = delta.*SI.e;

t_at = t.*(delta.*SI.e/SI.h_); 
dt_at = t_at(2)-t_at(1);

alpha = sigma*abs(cE_at).^2.*dt_at;


parfor i=1:length(I); 

    dW(i) = keldysh_rate_ln(lambda, 1, deltaSI, log(I(i)), m);
end;


W = zeros(size(I)); 
%sigma = 0;

for i=2:length(I); 
   W(i) = W(i-1) + dW(i)*(t(2)-t(1)) + alpha(i)*W(i-1);
end;

%v = cumtrapz(real(cE))*(t(2)-t(1)).*(delta/d).*SI.e/SI.m/m;




%jck = v(:).*W(:).*SI.e;


