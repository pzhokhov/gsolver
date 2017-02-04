function I = double_denom_integral(A,B,C,D);


T = B.*C - A.*D; 

L = log((A-B).*(C+D)./(A+B)./(C-D));

I = L./T; 
zT = abs(T)<1e-3; 
zC = abs(C)<1e-3; 
I(zT) = 2.*A(zT)./C(zT)./(B(zT).^2 - A(zT).^2);
I(zC) = 2./B(zC)./D(zC);