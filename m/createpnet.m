function net = createpnet(Xmin, Xmax, N, pw)

n0 = N/(1-oddroot(Xmax/Xmin,pw));
A = -Xmin/n0^pw;
net = A*((0:(N-1))-n0).^pw;


return;


function y = oddroot(x, p)
  if (mod(p+1,2)~=0) y = x.^(1/p); return; end;
  y = sign(x)*(abs(x)).^(1/p);
  return;
 
      
      