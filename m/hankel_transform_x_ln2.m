    function [x] = hankel_transform_x_ln2(N)


alpha = log(2);
x  = (1+exp(alpha))*exp(-alpha*(N-(0:(N-1))))/2;
