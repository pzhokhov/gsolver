function n = MeOH_n_omega(omega)
load const_SI;
lambda = 2*pi*SI.c./omega;
n = MeOH_n_lambda(lambda);
return;