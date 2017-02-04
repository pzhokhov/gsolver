
  
 Tmax = 2000; 
 dt = 0.1;
 t = (0:dt:Tmax-dt);

 tau = 300;  
 
 E0 = 1e-8; 
 E = E0.*exp(-((t-Tmax/2)/tau).^2).*cos(2*pi*(t-Tmax/2)./T0);


[d1,~,~,f1,g1,x] = H_atom_(sqrt(24.4/13.6), t, E); 
[d2,~,~,f2,g2,x] = H_atom_(2,               t, E);


alpha1 =  1.8076; 
alpha2 =  0.1924;
