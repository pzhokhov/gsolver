function [d, t, E, f_, g_, x] = H_atom_(Zeff, t, E, wb)
Nz = 256; 


if (nargin < 1) Zeff = 1; end;
if (nargin < 2);
       dt = 0.003; 
       Tmax = 400; 
       t = (0:dt:Tmax-dt);

       tau = 30;  
       T0 = 30;
       E0 = 10;
       E = E0.*exp(-((t-Tmax/2)/tau).^2).*cos(2*pi*(t-Tmax/2)./T0);
end; 

if (nargin < 4); wb = 1; end;

E = E/Zeff^3; 
t = t*Zeff^2; 
zmax = 20; 

Nt = length(t);

x  = (zmax/Nz : zmax/Nz : zmax)';
 

f = zeros(Nz, Nt); 
g = f; 

% f0 = exp(-z./2); 
% g0 = f0; 
% 
% f(:,1) = f0; 
% g(:,1) = g0; 
% 
f_ = f; g_ = g;
f_(:,1) = ones(Nz,1); 
g_(:,1) = ones(Nz,1); 


if (wb == 1); h = waitbar(0, 'Calculating again, please wait ...'); end;

dx = x(2)-x(1); 
%k1 = zeros(Nz,1); k2=k1; 
kappaf = zeros(Nz,1); gammaf = kappaf; gammag = kappaf; kappag = kappaf; 

%dt = 0; 

for nt=2:Nt; 
  dt = t(nt)-t(nt-1);
  kappaf(2) = 1; gammaf(2) = 0; 
  kappag(2) = 1; gammag(2) = 0; 
  for nz=2:(Nz-1);
      
    A = dt/dx^2 - dt*(1./x(nz) - 1)/2/dx;
    C = dt/dx^2 + dt*(1./x(nz) - 1)/2/dx;
    
    Bf = - 2*dt/dx^2 + dt*x(nz)*E(nt  )/2 + 1i; 
    Df = - 2*dt/dx^2 + dt*x(nz)*E(nt-1)/2 - 1i; 
    
    Ff = -A*f_(nz-1, nt-1) - C*f_(nz+1, nt-1) - Df*f_(nz, nt-1); 
    
    Gf = (A*kappaf(nz)+Bf); 
    kappaf(nz+1) = -C/Gf; 
    gammaf(nz+1) = (Ff-A*gammaf(nz))/Gf;
    
    
    Bg = - 2*dt/dx^2 - dt*x(nz)*E(nt  )/2 + 1i; 
    Dg = - 2*dt/dx^2 - dt*x(nz)*E(nt-1)/2 - 1i; 
    
    Fg = -A*g_(nz-1, nt-1) - C*g_(nz+1, nt-1) - Dg*g_(nz, nt-1); 
    
    Gg = (A*kappag(nz)+Bg); 
    kappag(nz+1) = -C/Gg; 
    gammag(nz+1) = (Fg-A*gammag(nz))/Gg;
  end; 
  f_(Nz, nt) = gammaf(Nz)/(1-kappaf(Nz)); 
  g_(Nz, nt) = gammag(Nz)/(1-kappag(Nz));   
  
  for nz=(Nz-1):-1:1;
    f_(nz, nt) = kappaf(nz+1)*f_(nz+1, nt) + gammaf(nz+1); 
    g_(nz, nt) = kappag(nz+1)*g_(nz+1, nt) + gammag(nz+1); 
  end; 
    
%     k1(2:Nz-1) = dt*2i*((f_(1:end-2, nt-1) - 2*f_(2:end-1,nt-1) + f_(3:end, nt-1))/dx^2 + (f_(3:end,nt-1)-f_(1:end-2,nt-1))/dx.*(1./x(2:end-1)-1) + f_(2:end-1,nt-1)/2*E(nt).*x(2:end-1));
%     k1(1)=k1(2); k1(Nz)=k1(Nz-1);
%     f1 = f_(:,nt-1) + k1; 
%     k2(2:Nz-1) = dt*2i*((f1(1:end-2)       - 2*f1(2:end-1)      + f1(3:end))/dx^2       + (f1(3:end)     -f1(1:end-2))/dx     .*(1./x(2:end-1)-1) + f1(2:end-1     )/2*E(nt).*x(2:end-1));
%     k2(1) = k2(2); k2(Nz)=k2(Nz-1);  
%     
%     f_(:,nt) = f_(:,nt-1) + 0.5*(k1+k2); 
%     
    if (wb == 1) h = waitbar(nt/Nt, h); end;
end; 

if (wb == 1); close(h); end;

for i=0:2; 
   %If(:,i+1) = squeeze(trapz(x, abs(f_).^2.*((x.^(i).*exp(-x))*ones(1,Nt))));
   %Ig(:,i+1) = squeeze(trapz(x, abs(g_).^2.*((x.^(i).*exp(-x))*ones(1,Nt))));
   If(:,i+1) = squeeze(sum(abs(f_).^2.*((x.^(i).*exp(-x))*ones(1,Nt))))*dx;
   Ig(:,i+1) = squeeze(sum(abs(g_).^2.*((x.^(i).*exp(-x))*ones(1,Nt))))*dx;
end;   

d = 0.5/Zeff*(If(:,3).*Ig(:,1) - Ig(:,3).*If(:,1))./(If(:,2).*Ig(:,1)+If(:,1).*Ig(:,2));

return; 
