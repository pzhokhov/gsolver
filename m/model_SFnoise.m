load const_SI light_v vacuum_permittivity
wl0 = 355e-9;
omega0 = 2*pi*light_v/wl0;
Omega = (-0.05:0.0005:0.0495);  %normalized fluctuation frequency
N     = length(Omega);
kt    = (0:0.5:100);        %normalized fluctuation transverse wavenumber  
a0    = 1e-3;
tau_  = [100 50 10]*1e-15;
for nt = 1:length(tau_);
tau   = tau_(nt);
    
beta  = 1e3;              %power to critical power ratio
f = 1e-3*[4:60];
for nf = 1:length(f);
kappa = 1./f(nf);              %initial beam curvature (normalized to inverse diffraction length
zR = pi*quartz_n_omega(omega0)*a0^2/wl0;
k = quartz_n_omega(omega0*(1+Omega)).*omega0.*(1+Omega)/light_v; 
k0 = quartz_n_omega(omega0)*omega0/light_v;
U0 = omega0*(Omega(N/2+1)-Omega(N/2-1))/(k(N/2+1)-k(N/2-1));
U = omega0*diff(Omega)./diff(k); U(end+1)=U(end);
d = k-omega0*Omega/U0-k0;  %unnormalized HO dispersion
T = U0*tau/zR;
zgr = abs(T./(U./U0-1));

[KT,OMEGA] = meshgrid(kt, Omega);
[KT,D]     = meshgrid(kt,d*zR); 
[KT,ZGR]   = meshgrid(kt, zgr);
E = ones(size(KT));

if (beta > 1)
 znlf = (sqrt(beta-1)-kappa)/(beta-1-kappa^2);    %!!!Achtung! This value of zmax iz valid only for beta > 1. Weak beam analysis should be performed with different value. 
else
 znlf = 1;
end
zmax = znlf;
Nz = 100;
z = 0:zmax/(Nz-1):0.99*zmax;
%z = [0, 0.03*zmax, 0.06*zmax, 0.1*zmax];
a2 = 1-2*kappa*z+(1-beta+kappa^2)*z.^2;
LAMBDA = zeros(size(KT,1),size(KT,2), length(z));
G      = zeros(size(LAMBDA));

for nz = 2:length(z);
Z = ones(size(KT))*z(nz); E_ = ones(size(KT)); E_(Z>ZGR)=0;
Theta = beta/a2(nz)*E_;
LAMBDA(:,:,nz) = sqrt((KT.^2./4./(1+OMEGA) - D - beta/a2(nz).*OMEGA - 1/a2(nz)).*(-KT.^2./4./(1+OMEGA) + D + beta/a2(nz).*(2+3.*OMEGA) + 1./a2(nz)));    
G(:,:,nz) = G(:,:,nz-1)+(z(nz)-z(nz-1)).*sqrt((KT.^2./4./(1+OMEGA) - D - Theta.*OMEGA - 1/a2(nz)).*(-KT.^2./4./(1+OMEGA) + D + Theta.*(2+3.*OMEGA) + 1/a2(nz)));
end;
maxg = squeeze(max(max(real(G)))); 
if (max(maxg(1:end-1)==maxg(2:end))==1)
    warning('maxg is ambigous');
end;
zt(nf,nt) = interp1(maxg, z, log(100))/znlf;  
end;
end;




