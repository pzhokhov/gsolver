function psi= qr_cohstate(phN, af, as, ap);

psi = zeros(phN^3,1);

for nf=0:phN-1; for ns=0:phN-1; for np=0:phN-1;
  n = 1 + nf+ ns*phN + np*phN^2;
 
  psi(n) = af^nf * as^ns * ap^np / sqrt(factorial(nf)*factorial(ns)*factorial(np));
             
end;end;end;

psi = psi/sqrt(psi'*psi); 