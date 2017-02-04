function I = bessel_exponent_integral(x,a)

I = zeros(size(x)); 
for na = 1:length(x(:)); 

 I(na) = integral(@(t)(exp(1i*a*t).*besselj(0,t)), 0, x(na));  
end;
