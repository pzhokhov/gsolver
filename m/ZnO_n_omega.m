function n = ZnO_n_omega(omega)

const_SI;
lambda = 2*pi.*(SI.c)./omega;
n = ZnO_n_lambda(lambda);