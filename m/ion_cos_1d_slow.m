function [W,W_] = ion_cos_1d_slow(a0, T, cep, alpha)

if nargin <= 2
 cep = 0;
end;

if nargin <= 3
  alpha = 1; 
end;

const_SI; 

w = 0.2123;

Nt = 8192; 
t_cyc = -7*T + (0:(Nt-1))*14*T/Nt;  t_cyc=t_cyc(:);
t = t_cyc./w*2*pi; 


dt = t(2)-t(1);

E = a0.*w.*real(exp(1i*t_cyc*2*pi - (t_cyc/T).^2+1i*cep));  %electric field

A = cumsum(E)*dt; %vector potential


cA = cumsum(cos(A))*dt*alpha;
sA = cumsum(sin(A))*dt*alpha;


x = (0:8191)*alpha.*dt.'/(8192./Nt);
nu =(1+alpha)/alpha.*(0.8:0.001:1.1);


[Nu,X] = meshgrid(nu,x);
%ba1_table = bessel_exponent_integral (phigrid, (1+alpha)/alpha); 
%ba2_table = besselx_exponent_integral(phigrid, (1+alpha)/alpha); 

ba1_table = zeros(length(X), length(nu));
ba2_table = zeros(length(X), length(nu));

for nnu = 1:length(nu)
    
  ba1_table(:,nnu) = cumtrapz(exp(1i*nu(nnu).*x).*besselj(0,x))   .*(x(2)-x(1)); 
  ba2_table(:,nnu) = cumtrapz(exp(1i*nu(nnu).*x).*besselj(0,x).*x).*(x(2)-x(1)); 

end


%E = (0:(Nt-1)).';

jpa = zeros(Nt,1); 
Ijc_c = zeros(Nt,1); 
Ijc_s = zeros(Nt,1); 

I = zeros(Nt,1);
I_ = I;
I__ = I;

tic_gp;




for nt2 = 2:Nt; 
    C = cA(nt2) - cA(1:nt2);
    S = sA(nt2) - sA(1:nt2);
    
    phiA = sqrt(C.^2 + S.^2); 
    
    %phiA = cumsum(ones(nt2,1))*dt*alpha; phiA = phiA(nt2) - phiA;
    
    beta = [alpha; -diff(phiA)/dt]; 
    
        
    P0 = besselj(0, phiA);
    %P1 = besselj(1, phiA); 

    
    %argphi = atan2(S,C); 
    
    tau = t(nt2)-t(1:nt2);  
    
    K = exp(1i*(1+alpha)*tau);
    M = real(K.*P0);
    
     
    
    %ba1 = interp2(X, Nu  ba1_table, phiA, (1+alpha)./beta); 
    ba1 = interp2(Nu, X, ba1_table, (1+alpha)./beta, phiA); 
    ba2 = interp2(Nu, X, ba2_table, (1+alpha)./beta, phiA); 
    
    
    dba1 = [0; -diff(ba1)];
    dba2 = [0; -diff(ba2)];
    
%     ba2 = zeros(nt2,1);
%     ba2(nt2) = interp1(phigrid, ba2_table, phiA_(nt2));
%     ba2(nt2-1) = 0; 
%     ba2(1:(nt2-2)) = -conj(interp1(phigrid, ba2_table, -phiA_(1:(nt2-2))));
%     
%    dba2 = [0; diff(ba2)]/alpha/dt;
        
    dE =  [0; -diff(E(1:nt2))];
    
    M1_ = real(exp(1i*(1+alpha)*tau - 1i*phiA.*(1+alpha)./beta).*dba1)./beta;
    M2_ = real(exp(1i*(1+alpha)*tau - 1i*phiA.*(1+alpha)./beta).*(dba2 - beta.*tau.*dba1)./beta.^2);
    
    I(nt2)  =  2*trapz(M(1:nt2).*E(1:nt2)).*dt.^2;
    I__(nt2) = 2*sum((M1_.*E(1:nt2).*dt + M2_.*dE));
    I_(nt2) =  2*sum(M1_(1:nt2).*E(1:nt2));
    toc_gp(nt2/Nt);
end

 
 
W   = cumtrapz(I.*E);
W_  = cumtrapz(I_.*E);
W__ = cumtrapz(I__.*E);
%jc = (2*pi)./2.*delta.*dt.^2*(cumtrapz(E.*Ijc_s(:))); %.*cos(A) - cumtrapz(E.*Ijc_c(:)).*sin(A));

%jpa = I;
% t = w*t; 

%W = -cumtrapz(E.^2.*(cos(A)-1)).*alpha.^2.*dt + E.^2.*(1+alpha);
return;

