function k2 = quartz_GVD_omega(omega)

const_SI;

delta = 1e-4;
omega_ = omega(:);
o = omega(:)*[1-delta, 1, 1+delta];
n = quartz_n_omega(o); 
k = o.*n./SI.c;

k2 = (k(:,1)+k(:,3)-2*k(:,2))./(delta.*omega_).^2;
k2 = reshape(k2, size(omega));