function psi= qr3_cohstate(phN, af, as, ap, aas)

psi = zeros(phN^4,1);

for nf=0:phN-1; for ns=0:phN-1; for np=0:phN-1; for nas=0:phN-1;
  n = 1 + nf+ ns*phN + np*phN^2 + nas*phN^3;
 
  psi(n) = af^nf * as^ns * ap^np * aas^nas / sqrt(factorial(nf)*factorial(ns)*factorial(np)*factorial(nas));
             
end;end;end;end;

psi = psi/sqrt(psi'*psi); 