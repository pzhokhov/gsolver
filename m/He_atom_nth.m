function [khi, khi1, khi2] = He_atom_nth(lambda);
 load const_SI; 
 a0 = 5.3e-11;
 t0 = 2.4189e-17;
 Eat = 5.1421e11;
 h = waitbar(0, 'Wait...'); 
 for nl = 1:length(lambda); 
 
  T0 =  lambda(nl)./SI.c./t0; 
  
  Tmax = 2000; 
  dt = 0.1;
  t = (0:dt:Tmax-dt);

  tau = 300;  
 
  E0 = 1e-8; 
  E = E0.*exp(-((t-Tmax/2)/tau).^2).*cos(2*pi*(t-Tmax/2)./T0);


 [d1,~,~,f1,g1,x] = H_atom_(sqrt(24.4/13.6), t, E, 0); 
 [d2,~,~,f2,g2,x] = H_atom_(2,               t, E, 0);


 alpha1 =  1.81; 
 alpha2 =  0.19;

 Psi1 = SI.Loschmidt_number*SI.e*a0*d1;
 Psi2 = SI.Loschmidt_number*SI.e*a0*d2;
 
 khi1(nl) = max(abs(fft(d1)))./max(abs(fft(E)))*SI.Loschmidt_number*SI.e*a0/SI.epsilon0./Eat;
 khi2(nl) = max(abs(fft(d2)))./max(abs(fft(E)))*SI.Loschmidt_number*SI.e*a0/SI.epsilon0./Eat;

 khi(nl) = alpha1*khi1(nl) + alpha2*khi2(nl); 
 waitbar(nl/length(lambda), h); 
 end;