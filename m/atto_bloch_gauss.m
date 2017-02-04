function [d,rho,t,E] = atto_bloch_gauss(cep, F)

d13 = 0;
Nt = 4096; 
tmin = -100;
tmax = 100;
 
 t = tmin + (tmax-tmin)*(1:Nt)/Nt;
 T = 6;
 const_SI;
 w0 = 2*pi*SI.c/800e-9/1e15; 
 
 E = F*real(exp(-(t/T).^2+1i*t+1i*  cep));
%  
%  load('attotrans_data.mat'); 
%  %t = attotrans_data(:,1)*w0; 
%  E = attotrans_data(:,2); 
%  %Nt = size(attotrans_data,1);
%  E = interp1(attotrans_data(:,1)*w0, attotrans_data(:,2), t, 'cubic', 0); 
%  E = E./max(abs(E(:)))*F;  %E(isnan(E))=0;
%  E = E.*exp(-(t/max(attotrans_data(:,1))/w0/0.9).^8);
%  fE = fft(E); 
%  fE(1:(length(fE)/2))            = fE(1:(length(fE)/2))           *exp( 1i*cep); 
%  fE((length(fE)/2+1):length(fE)) = fE((length(fE)/2+1):length(fE))*exp(-1i*cep); 
% % fE(abs(fE)<1e-3*max(abs(fE)))=0;
% 
%  E = ifft(fE);
 w1 = 6.47;
 w2 = 8.3882;
 Q2 = 100; 
 Q1 = 60;
 
%  bloch_sys.N = 3;
%  bloch_sys.en = [0 w1 w2];
%  bloch_sys.d  = [0 1 d13; 1 0 0; d13 0 0];
%  bloch_sys.G2 = [0 w1/Q2 w2/Q2; w1/Q2 0 0; w2/Q2 0 0];
%  bloch_sys.G1 = [0 0 0;0 0 0;0 0 0];
%  
 bloch_sys.N = 2;
 bloch_sys.en = [0 w1];
 bloch_sys.d  = [0 1;1 0 ];
 bloch_sys.G2 = [0 w1/Q2;w1/Q2 0];
 bloch_sys.G1 = [0 0;0 0];
 


 rho0 = zeros(bloch_sys.N,bloch_sys.N); rho0(1,1)=1;  %zero temperature; 
 yin = rho0(:); 
 
 [~, yout] = ode45(@(x,y)bloch_rhs(x,y,bloch_sys,t,E), t, yin); 
 
 rho = reshape(yout, Nt, bloch_sys.N, bloch_sys.N);
 
 for i=1:Nt; 
  d(i) = trace(squeeze(rho(i,:,:))*bloch_sys.d); 
 end;

 %d = cat(1,d(:), zeros(Nt,1)); 
 return; 
 

function [R]=bloch_rhs(x, y, p, t, E)
E_ = interp1(t,E,x);

rho  = reshape(y, p.N, p.N); 
H = diag(p.en) - (p.d*E_); 
drho = -1i*(H*rho - rho*H) - rho.*p.G2;

for i=1:p.N
 for k=1:p.N
     drho(i,i) = drho(i,i) + (p.G1(i,k)*rho(k,k) - p.G1(k,i)*rho(i,i));
 end;
end;

R = drho(:);
return; 