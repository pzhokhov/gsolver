function [Eout] = atto_bloch_prop(cep, F,T)

d13 = 0;
Nt = 2048; 
tmin = -100;
tmax = 100;
 
 t = (tmin + (tmax-tmin)*(1:Nt)/Nt).';
% T = 6;
 const_SI;
 w0 = 2*pi*SI.c/800e-9/1e15; 
 
 z = [0:0.01:0.5]; 
 
 E0 = F*(exp(-(t/T).^2+1i*t+1i*  cep));
 %h = waitbar(0, 'Calculating, please wait...');
 
 [~, Eout] = ode113(@(x,y)bloch_prop_RHS(x,t,y,z), z, E0); 
 
 Eout = Eout.'; 
 
% close(h);
 return; 
 
 
 function [dE] = bloch_prop_RHS(z, t, E, zspan);
%  E = ifft(fE);


 Nt = length(t); 
 w1 = 6.47;
 w2 = 8.3882;
 Q2 = 30; 
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
 
 [~, yout] = ode45(@(x,y)bloch_rhs(x,y,bloch_sys,t,real(E)), t, yin); 
 
 rho = reshape(yout, Nt, bloch_sys.N, bloch_sys.N);
 
 d=zeros(Nt,1);
 for i=1:Nt; 
  d(i) = trace(squeeze(rho(i,:,:))*bloch_sys.d); 
 end;

 d = d+0.01.*E.^3; 

% d = 0.1.*E.^3; 
 
 w  = ifftshift(cfreq(t)).'; 
 fd = fft(d); 
  
 fdE = -1i.*w.*fd; 
 fdE(w<0) = 0; 
 fdE(w>0.7*max(w))=0;
 dE = ifft(fdE);
 dE = dE.*exp(-(t(:)/0.8/max(t)).^10); 
 
 %waitbar((z-min(zspan))/(max(zspan)-min(zspan)),h);
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
