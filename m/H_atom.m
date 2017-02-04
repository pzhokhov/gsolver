Nz = 256; 

Tmax = 200; 

zmax = 30; 

dt =0.01; 
tau = 20;  
T0 = 5;
t = (0:dt:Tmax-dt);
Nt = length(t);
E0 = 0.0000; 

E = E0.*exp(-((t-Tmax/2)/tau).^2).*cos(2*pi*(t-Tmax/2)./T0);

[x,V] = Hankel_transform_bessel_x(Nz);
z  = x*zmax;
kz = x*2*pi*V./zmax; 

 

f = zeros(Nz, Nt); 
g = f; 

f0 = exp(-z./2); 
g0 = f0; 

f(:,1) = f0; 
g(:,1) = g0; 


h = waitbar(0, 'Calculating again, please wait ...');



load c.mat;
c = c(1,1:Nz+1);

[Jn,Jm] = meshgrid(c(1:Nz),c(1:Nz));
C = (2/c(Nz+1))*besselj(0,Jn.*Jm/c(Nz+1))./(abs(besselj(1,Jn)).*abs(besselj(1,Jm)));
m1 = (abs(besselj(1,c(1:Nz))))';

for nt=2:Nt; 
    
    hf = ((C*(f(:,nt-1)./m1)));
    %hf_ = Hankel_transform_bessel(f(:,nt-1),1);
    
    hf = hf .* exp(-2i.*kz.^2.*dt); 
    hf = (C*(hf)).*m1;
    vf = 1./z + z.*E(nt) - 0.48;
    
    f(:,nt) = hf.*exp(1i*dt.*vf); 
    
    hg = ((C*(g(:,nt-1)./m1)));
    %hf_ = Hankel_transform_bessel(f(:,nt-1),1);
    
    hg = hg .* exp(-2i.*kz.^2.*dt); 
    hg = (C*(hg)).*m1;
    vg = 1./z - z.*E(nt);
    
    g(:,nt) = hg.*exp(1i*dt.*vg); 
   
    h = waitbar(nt/Nt, h); 
end; 

close(h); 