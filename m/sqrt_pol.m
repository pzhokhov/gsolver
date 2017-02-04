function f = sqrt_pol(x,b);

N = 20      ; 

f = ones(size(x)) + x/2;
if (nargin == 1) 
    b = 0;
end;
if (b == 1)
    T = 1/2;
end;
for i=2:N;
    M  = fac2(2*i-3);
    t = -(-1/2)^i*M/factorial(i);
    f = f+t*x.^i;
    if (b==1)
        T = [T;t];
    end;
end;

if (b==1)
  T
end;

  

function d = fac2(N);
d = prod(N:-2:1);