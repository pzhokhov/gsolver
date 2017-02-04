function S = S_integral_1d(a,b1,b2)

if (size(a)~=size(b1))
    error('Size of a is not equal to the size of b1!');
end;
if (size(a)~=size(b2))
    error('Size of a is not equal to the size of b2!');
end;

S = zeros(size(a)); 


pmaxx = 2; 
px = -pmaxx + 2*pmaxx*(0:511)./512;


S = zeros(size(a));
parfor nx = 1:length(px);
    K = exp(-1i.*px(nx).^2.*a)./(1+(px(nx)).^2).^2; 
    S = S + K; 
end;
S = S.*(px(2)-px(1));



