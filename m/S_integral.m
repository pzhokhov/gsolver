function S = S_integral(a,b); 

if (size(a)~=size(b))
    error('Size of a is not equal to the size of b!');
end;

a = abs(a); 
b = abs(b); 
a_ = a(:);
b_ = b(:); 


S_ = zeros(1,length(a_)); 




%S = pi^(3/2).*sqrt(abs(a)).*exp(-1i*a-3*pi/4) + pi^2.*(0.5+1i.*abs(a)).*(1-erfz(exp(1i.*pi/4).*sqrt(abs(a))));
%return; 

%S = i./2.*(1i*pi./a).^(3/2).*exp(-1i*a); return;

a_ = abs(a_); 

parfor n = 1:length(a_); 
%     hx = min(0.01, a(n)/100); 
%     minx = min(b(n)-5*a(n), -10);
%     maxx = max(b(n)+5*a(n), +10); 
%     

       if (b_(n) <= 1e-5)
          if (a_(n) <=1e-5);
            S_(n) = 2*pi^2; continue;
          else
            hx = 1e-3*(1 + 1./a(n)).^(-1);
            minx = -50; maxx = 50;
            x = minx : hx : maxx;
    
            T = 2.*x.^2./(1 + (x./a(n)).^2).^2.*(exp(-1i*x.^2)); 
            S_(n)  = 2*pi*sum(T)*hx./(a_(n)).^(3).*exp(-1i*a(n).^2);   
          end;
       else
         if (a(n) <= 1e-5)
             S_(n)=0;
         else
            hx = 0.01*(1 + 1./a(n)).^(-1);
            minx = -50; maxx = 50;
            x = minx : hx : maxx;
    
            T = exp(-1i*x.^2).*x./(1 + ((x-b_(n))./a(n)).^2); 
            S_(n)  = 2*pi*sum(T)*hx./b_(n)./(a_(n)).*exp(-1i*a(n).^2);
         end;
       end; 
end;


% parfor n = 1:numel(a);  
%     if (a(n) < sqrt(pi))
%         if (a(n)==0 && b(n) < 1e-5) 
%             S_(n) = 2*pi.^2;
%             continue; 
%         end;
%         if (a(n)==0 && b(n) > 1e-5)
%             S_(n) = 0;
%             continue;
%         end; 
%                 
%         if (b(n)/a(n) <= 1e-5); 
%             small a, b=0 case
%             hx = 0.01;
%             minx = -200; maxx = 200;
%             x = minx : hx : maxx;
%             
%             T = 2.*x.^2./(1 + x.^2).^2.*exp(-1i*a(n).^2.*x.^2); 
%             S_(n)  = 2*pi*sum(T).*hx.*exp(-1i*a(n).^2);   
%             continue;
%         else
%             hx = 0.1;
%             minx = -1000+b(n)/a(n); maxx = 1000+b(n)/a(n);
%             x = minx : hx : maxx;
%             
%             T = exp(-1i*a(n)^2.*x.^2).*x./(1 + (x-b(n)./a(n)).^2); 
%             S_(n)  = 2*pi*sum(T)*hx./b_(n).*exp(-1i*a(n).^2).*a(n);
%             continue; 
%         end;
%     else
%         if (b(n) <= 1e-5)
%             minx = -100; maxx = 100; 
%             hx = 0.01;
%             x = minx : hx : maxx;
%             T = 2.*x.^2./(1 + (x./a(n)).^2).^2.*exp(-1i*x.^2); S_(n)  = sum(T)*hx./b_(n).*exp(-1i*a(n).^2).*a(n);
%             S_(n)  = 2*pi*sum(T)*hx./(a_(n)).^(3).*exp(-1i*a(n).^2);   
%         else
%             minx = -100; maxx = 100; 
%             hx = 0.01;
%             x = minx : hx : maxx;
%             T = exp(-1i*x.^2).*x./(1 + ((x-b_(n))./a(n)).^2); 
%             S_(n)  = 2*pi*sum(T)*hx./b_(n)./(a_(n)).*exp(-1i*a(n).^2);
%         end;
%     end;
% end;


S = reshape(S_, size(a)); 

