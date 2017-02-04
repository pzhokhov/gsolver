function n = ZGP_no_lambda(lambda)

lambda = lambda/1e-6;

n2 = 4.61511 + 5.12798*lambda.^2./(lambda.^2 - 0.13624) + 2.16936*lambda.^2./(lambda.^2-900);

n = n2.^(0.5);