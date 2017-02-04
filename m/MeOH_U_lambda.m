function k2 = MeOH_U_lambda(lambda);

load const_SI ;
omega = 2*pi*SI.c./lambda;

k2 = MeOH_U_omega(omega);