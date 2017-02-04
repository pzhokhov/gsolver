function B = besselsum4_approx(X1, Y1, X2, Y2, X3, Y3, X4, Y4, nmax); 

if (nargin == 8)
    nmax = 8; 
end; 
    x =  sqrt(X1.^2 + Y1.^2); phi1 = atan2(Y1, X1);
    
    
    phi = 5*pi/4; A = -0.0454; f = 4;
   
    B = zeros(size(x));
%     for i=1:length(x(:))
%      if (x(i)<0.2) B(i) = besselj(0,x(i)).^4; 
%      else B(i) = 1./(1 + (x(i)/0.395).^2).^0.425; %+A*(cos(4*x(i) - phi))./(x(i).^(5/2)) + 
%      end;
%     end;
    B = besselj(0,4*x);
end