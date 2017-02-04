function S = S_integral_3d(a,b1,b2)

if (size(a)~=size(b1))
    error('Size of a is not equal to the size of b1!');
end;
if (size(a)~=size(b2))
    error('Size of a is not equal to the size of b2!');
end;

S = zeros(size(a)); 


pmaxx = 2; 
%px = -pmaxx + 2*pmaxx*(0:511)./512;
px = pmaxx*(0:511)./512;


S = zeros(size(a));
St = zeros(length(a), length(px));
parfor nx = 1:length(px);
      p = px(nx);
%     K = exp(-1i.*p^2.*a).*p./((b1-b2).*(1+p.^2-b1.*b2)).*log((1+(p+b1).^2)./(1+(p+b2).^2));
%     
%     K(abs(b1-b2)<1e-5) = exp(-1i.*p^2.*a(abs(b1-b2)<1e-5)).*p./(1+(p-b1(abs(b1-b2)<1e-5)).^2)./b1(abs(b1-b2)<1e-5); 
%     K(abs(b1-b2)<1e-5 & (b1 < 1e-5)) = 2.*exp(-1i.*p^2.*a((abs(b1-b2)<1e-5) & (b1 < 1e-5 ))).*p.^2./(1+p.^2).^2;
%     
%    K = exp(-1i.*p^2.*a).*p.^2./(1+(p+b1).^2)./(1+(p+b2).^2);
       
    D = double_denom_integral(2*p*b1, 1+p^2+b1.^2, 2*p*b2, 1+p^2+b2.^2); 
    K = exp(-1i.*p^2.*a).*p.^2.*D;
    
    S = S + K;   
end;
S = S.*(px(2)-px(1))*pi*2;



