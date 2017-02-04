function k2 = quartz_GVD_lambda(lambda);

const_SI;
omega = 2*pi*SI.c./lambda;

k2 = quartz_GVD_omega(omega);