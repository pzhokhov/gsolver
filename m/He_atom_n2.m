function [khi, khi1, khi2, d1_, d2_] = He_atom_n2(I);
 load const_SI; 
 a0 = 5.3e-11;
 t0 = 2.4189e-17;
 Eat = 5.1421e11;
 h = waitbar(0, 'Wait...'); 
 
 lambda = 800e-9; 
 C = sqrt(2*sqrt(SI.mu0/SI.epsilon0));
 
 for nl = 1:length(I); 
 
  T0 =  lambda./SI.c./t0; 
  
  Tmax = 700; 
  dt = 0.01;
  t = (0:dt:Tmax-dt);

  tau = 100;  
 
  E0 = sqrt(I(nl))/C/Eat;
  E = E0.*exp(-((t-Tmax/2)/tau).^2).*cos(2*pi*(t-Tmax/2)./T0);


 [d1,~,~,f1,g1,x] = H_atom_(1, t, E, 1); 
 [d2,~,~,f2,g2,x] = H_atom_(1,               t, E, 1);


 alpha1 =  1.81; 
 alpha2 =  0.19;

 Psi1 = SI.Loschmidt_number*SI.e*a0*d1;
 Psi2 = SI.Loschmidt_number*SI.e*a0*d2;
 
 khi1(nl) = max(abs(fft(d1)))./max(abs(fft(E)))*SI.Loschmidt_number*SI.e*a0/SI.epsilon0./Eat;
 khi2(nl) = max(abs(fft(d2)))./max(abs(fft(E)))*SI.Loschmidt_number*SI.e*a0/SI.epsilon0./Eat;

 khi(nl) = alpha1*khi1(nl) + alpha2*khi2(nl); 
 
 d1_(:,nl) = d1; 
 d2_(:,nl) = d2; 
 
 waitbar(nl/length(I), h); 
 end;
 
 close(h);