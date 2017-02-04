function B = besselsum4_interp(X1, Y1, X2, Y2, X3, Y3, X4, Y4, nmax); 

global besselsum4_bs;
global besselsum4_x; 

if (isempty(besselsum4_x))
    load besselsum4;
    besselsum4_bs = bs4; 
    besselsum4_x = x; 
end;

bs4 = besselsum4_bs; 
x   = besselsum4_x; 

if (nargin == 8)
    nmax = 8; 
end;
    load besselsum4; 
    rho1 =  sqrt(X1.^2 + Y1.^2); phi1 = atan2(Y1, X1);
    
    B = interp1(x,bs4, rho1, 'linear', 0); 
end