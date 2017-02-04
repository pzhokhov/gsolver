function [d, E, t, rho] = atto_bloch(cep, I, e12, d12, Q2, Efun)



    
const_SI;

 peaks = cep; 
 Ncep = size(peaks,2); 
 
 if (nargin <= 5)
     Efun = @(a,b)atto_E_peaks(a,b,8); 
 end;
 
 [E,t] = Efun(peaks, I);

%w0 = 2*pi*SI.c/800e-9/1e15;
w12 = e12.*SI.e/SI.atomic_energy;


%d12 = (3*SI.h_.*f12/2/SI.m./w12/1e15).^(0.5)*SI.e;
Emax = sqrt(2*I*sqrt(SI.mu0/SI.epsilon0))/SI.atomic_field;
E = E./max(abs(E)).*Emax;

Nlev = numel(w12); 
Nt = length(t); 
Ncep = size(cep, 2);
w12 = reshape(w12, 1, Nlev); 
d12 = reshape(d12, 1, Nlev); 

%Omega = d12*Emax;

%gamma12=w12./Q2; 

 
 
 t = t*1e-15/SI.atomic_time; %convert t into the atomic units of time

 
 
 
 
 bloch_sys.N = Nlev+1;
 bloch_sys.en = [0 w12];
 
 bloch_sys.d = zeros(bloch_sys.N, bloch_sys.N);
 bloch_sys.G2 = zeros(bloch_sys.N, bloch_sys.N);
 bloch_sys.G1 = zeros(bloch_sys.N, bloch_sys.N);
 
 for i=2:Nlev+1; 
       bloch_sys.d(1,i)=d12(i-1); 
       bloch_sys.d(i,1)=d12(i-1); 
       bloch_sys.G2(1,i)=w12(i-1)/Q2(i-1); 
       bloch_sys.G2(i,1)=w12(i-1)/Q2(i-1); 
 end;
       
        
%  bloch_sys.d  = [0 1 d13; 1 0 0; d13 0 0];
%  bloch_sys.G2 = [0 w1/Q2 w2/Q2; w1/Q2 0 0; w2/Q2 0 0];
%  bloch_sys.G1 = [0 0 0;0 0 0;0 0 0];
%  
%  bloch_sys.N = 2;
%  bloch_sys.en = [0 w1];
%  bloch_sys.d  = [0 1;1 0 ];
%  bloch_sys.G2 = [0 w1/Q2;w1/Q2 0];
%  bloch_sys.G1 = [0 0;0 0];
 

for ncep = 1:Ncep;   
 
 dc = zeros(Nt,1);
    
 %cE = imag(hilbert(E0)*exp(1i*(pi*0.5+cep(ncep))));
 cE = E(:,ncep); 
 
 rho0 = zeros(bloch_sys.N,bloch_sys.N); rho0(1,1)=1;  %zero temperature; 
 yin = rho0(:); 
 
% [~, yout] = ode45(@(x,y)bloch_rhs(x,y,bloch_sys,t,cE), t, yin); 
 
% rho = reshape(yout, Nt, bloch_sys.N, bloch_sys.N);
 rho = zeros(Nlev+1, Nlev+1, Nt); 
 rho(:,:,1) = rho0;
 dt = t(2)-t(1); 
 tic_gp; 

 for nt=2:Nt;
    H = diag(bloch_sys.en) - (bloch_sys.d*E(nt)); 
    U = expm(1i*dt*H); 
    Ud = expm(-dt*bloch_sys.G2); 
    rho(:,:,nt) = U*rho(:,:,nt-1)*U'; 
    rho(:,:,nt) = rho(:,:,nt)*Ud;
 
    dc(nt) = trace(squeeze(rho(:,:,nt))*bloch_sys.d); 
    toc_gp(nt/Nt);
 end;

 d(:,ncep)=dc;
 
end;

 %d = cat(1,d(:), zeros(Nt,1)); 
 return; 
 

function [R]=bloch_rhs(x, y, p, t, E)
E_ = interp1(t,E,x);

rho  = reshape(y, p.N, p.N); 
H = diag(p.en) - (p.d*E_); 
drho = -1i*(H*rho - rho*H) - rho.*p.G2;

% for i=1:p.N
%  for k=1:p.N
%      drho(i,i) = drho(i,i) + (p.G1(i,k)*rho(k,k) - p.G1(k,i)*rho(i,i));
%  end;
% end;

R = drho(:);
return; 
