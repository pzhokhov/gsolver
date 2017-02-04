function [W,L,t] = ion_aniso_(gamma, T, delta)

const_SI; 
lambda = 800e-9; 

Nt = 2048; 
t = -40 + (0:(Nt-1))*80/Nt;  t=t(:);
t = t*1e-15/SI.atomic_time;
T = T*1e-15/SI.atomic_time;

dt = t(2)-t(1); 


w = 2*pi*SI.c/lambda/SI.atomic_frequency; 


mx = 0.9;
my = 0.6; 

delta = delta*w; 


A0 = sqrt(mx*delta)./gamma;

E = real(exp(1i*w*t - (t/T).^2));  %electric field
fE = fft(E);  fE(abs(fE)<1e-3*max(abs(fE)))=0; 
E = ifft(fE);

A = cumsum(E)*dt; %vector potential
E = E.*A0./max(A); E0 = max(E); 
A = A./max(A).*A0; 

% Ai = cumsum(A)*dt;
% A2i=cumsum(A.^2)*dt; 
% 
% E1E2 = E*E.';
% AI = Ai*ones(1,Nt); AI=-AI+AI.';
% A2I = A2i*ones(1,Nt); A2I=-A2I+A2I.'; 
% 
% Xi = t*ones(1,Nt); Xi=-Xi+Xi.';
% 
% % for nt=1:length(Nt);
% %  S(:,nt) = S_integral(delta*Xi(:,nt), sqrt(delta)*AI(:,nt)); 
% % end; 
% x1 = -5:0.02:5; x1=x1(:);
% x2 = -5:0.02:5; x2=x2(:).'; 
% 
% X1 = x1*ones(1,length(x2)); 
% X2 = ones(length(x1),1)*x2; 
% 
% %S0 = S_integral(X1, X2); 
% 
% 
% %S = interp2(X2, X1, S0, sqrt(delta)*(AI), delta*(Xi), 'linear', 0); 
% 
% %S = pi^2/2; S_integral(delta*Xi, zeros(size(Xi))); %sqrt(delta)*AI);
% 
% %load S_;
% 
% 
% sqXi = sqrt(abs(Xi)); 
% Ain = (abs(Xi.*delta));
% Bin = sqXi.*(ones(Nt,1)*A(:).') - AI./sqXi;
% Bin(Xi==0) = 0; 
% 
% %S = interp2(B_,A_,S_, abs(Bin), Ain, 'linear', 0); 
% 
% S = 2*(pi./Ain./1i).^(1/2)./(1 + (Bin./Ain).^2).^2.*exp(-1i.*Ain); 
% S(abs(Ain)<3*pi) = 0;
% 
% expA2I = exp(-1i.*(A2I - AI.^2./Xi)); expA2I(Xi==0)=1; 
% % hb = waitbar(0, 'Calculating S integral...'); 
% % for nt2 = 1:Nt; 
% %     parfor nt1 = 1:(nt2-1); 
% %         xi = Xi(nt1, nt2); 
% %         S(nt1, nt2) = S_integral(delta*xi, sqrt(xi).*A(nt2) - (Ai(nt2)-Ai(nt1)/xi));
% %         expA2I(nt1, nt2) = exp(-1i*(A2I(nt1, nt2)-AI(nt1, nt2)./xi)-i*xi*delta); 
% %     end;
% %     waitbar(nt2/Nt, hb); 
% % end;
% %close(hb); 
% % S = S_intergral(delta*Xi, sqrt(Xi).*(A.*ones(1,Nt) - AI
% % 
% % tA = exp(-sqrt(delta)./abs(A*ones(size(A)).')); tA(~(Xi==0))=0;
% % 
% % S = pi^2/2*exp(1i*delta*Xi)+tA;
% 
% 
% %expA2I = exp(-1i*(A2I-AI.^2./Xi)); expA2I(Xi==0)=1;  
% 
% I = E1E2.*expA2I.*S; 
% 



% 
% W = zeros(size(t)); W(1)=I(1,1); 
% 
% for nt=2:Nt; 
%     W(nt) = W(nt-1)+real(sum(I(nt,1:nt)))+real(sum(I(1:(nt-1),nt)));
% end; 
% 
% 
% W = W*dt;
% 
% 
% return; 

pmaxx = sqrt(mx*delta); 
px = -pmaxx + 2*pmaxx*(0:1023)./1024;

L = zeros(Nt, length(px));
parfor nx = 1:length(px);
    en = delta + (px(nx)+A).^2;
   
    K = exp(-1i*cumsum(en)*dt)./(delta + (px(nx) + A).^2); 
    L(:, nx) = dt*cumsum(K(:).*E(:)); 
end;
      

W = sum(abs(L).^2, 2)*dt; 
