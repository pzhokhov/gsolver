    function [x, alpha] = hankel_transform_x(N, al)
    
if (nargin == 1)
alpha = fsolve(@(y)(y+log(1-exp(-y))/(N-1)),1);
else
alpha = al;
end;

x0 = (1+exp(alpha))*exp(-alpha*N)/2;
x  = x0*exp(alpha*(0:(N-1))).';
