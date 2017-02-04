function [w] = ion_cos_1d_av(gamma, delta, alpha)


A0 = 1./gamma/sqrt(alpha);
w = 1;

E0 = A0./w; 

    S = 0;
    deltat = delta.*(1 + alpha.*A0.^2./4); 
    for k=1:100; 
            p = w*(2*k+1)-deltat;
            S = S+(p>0).*besselj(k, alpha*delta.*A0.^2/8/w).^2*sqrt(p);
    end; 
    w = sqrt(2)*pi*E0.^2./sqrt(alpha.*delta).*S;    
return;

