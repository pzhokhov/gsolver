function [W, E, t] = calculate_plasma_keldysh_a0_unitless(a0, T, al, w)



T = T*2*pi/sqrt(log(4)); 

Nt = 8192; 
t = -6*T/w + 12*T/w*(0:(Nt-1))/Nt;


Emax = w.*a0; 
Emax_at = a0*w;

cE =       Emax.* exp(-(w.*t/T).^2+ 1i*w.*t);
cE_at = Emax_at.* exp(-(w.*t/T).^2+ 1i*w.*t);


E = real(cE); 


parfor i=1:length(t); 
    dW(i) = keldysh_rate_unitless(abs(cE(i)), al, w);
end;


W = zeros(size(t)); 
%sigma = 0;

for i=2:length(t); 
   W(i) = W(i-1) + dW(i)*(t(2)-t(1));
end;

%v = cumtrapz(real(cE))*(t(2)-t(1)).*(delta/d).*SI.e/SI.m/m;




%jck = v(:).*W(:).*SI.e;


