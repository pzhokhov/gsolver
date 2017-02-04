function U = quartz_U_omega(omega)

const_SI;

delta = 1e-5;
omega_ = omega(:);
o = omega(:)*[1-delta, 1+delta];
n = quartz_n_omega(o); 
k = o.*n./SI.c;
U = ((k(:,2)-k(:,1))./(delta.*omega_)).^(-1);
U = reshape(U, size(omega));