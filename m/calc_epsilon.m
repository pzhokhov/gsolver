function E = calc_epsilon(filename, n, n0, n2, lambda0)

plasma_filename = [filename, '.plasma'];

I = abs(load_sectiondump(filename,n)).^2;
ro = load_sectiondump(plasma_filename, n);

const_SI;

w = 2*pi*light_v/lambda0;
roc = w^2*vacuum_permittivity*n0^2*e_mass/e^2;

E = n0^2*(1-ro/roc)+2*n0*n2*I - ro/roc;
