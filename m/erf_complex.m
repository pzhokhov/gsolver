function e = erf_complex(z)

z_ = z(:); 
e = zeros(size(z)); 

Nsum = 50;
for nz=1:length(z_);
    cz = z_(nz); 
    %if (imag(cz)==0) e_(nz) = erf(cz); 
    %else
        if (abs(cz)<3)
            S = 0;
            M = 1;
            for ns = 1:Nsum;
                S = S + cz/(2*ns-1)*M;
                M = M*(cz^2)*(-1)/(ns);
%                S= S+(-1)^(ns-1)*cz^(2*ns-1)/factorial(ns-1)/(2*ns-1);
            end;
            e(nz) = S*2/sqrt(pi);
        else
            S = 1;
            for ns = 2:Nsum;
                S = S + (-1)^dfac(2*ns-3)/(2*cz^2)^(ns-1);
            end;
            e(nz) = 1- S*exp(-cz^2)/cz/sqrt(pi);
        end;
    %end;
end;



function f = dfac(n)

f=zeros(size(n)); 
for i=1:numel(n);
 if (mod(n(i),2)==0)
  f(i) = prod(2:2:n(i))
 else
  f(i) = prod(1:2:n(i));
 end;
end; 
