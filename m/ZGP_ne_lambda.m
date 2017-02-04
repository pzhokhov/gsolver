function n = ZGP_ne_lambda(lambda, theta)

lambda = lambda/1e-6;
if (nargin < 2)
    theta = pi/2; 
end; 

n2 = 4.69874 + 5.27924*lambda.^2./(lambda.^2 - 0.14339) + 2.09861*lambda.^2./(lambda.^2-900);

ne = n2.^(0.5);
no = ZGP_no_lambda(lambda*1e-6);

n = 1./sqrt(sin(theta).^2./ne.^2 + cos(theta).^2./no.^2); 

