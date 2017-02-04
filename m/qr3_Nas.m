function A = qr3_Nas(phN)

A = zeros(phN^4);

for nf1=0:phN-1; for nf2=0:phN-1; for ns1=0:phN-1; for ns2=0:phN-1; for np1=0:phN-1; for np2=0:phN-1; for nas1=0:phN-1; for nas2=0:phN-1;
  n1 = 1 + nf1 + ns1*phN + np1*phN^2 + nas1*phN^3;
  n2 = 1 + nf2 + ns2*phN + np2*phN^2 + nas2*phN^3;
  
  
  A(n1,n2) = nas1*not(ns1-ns2)*not(np1-np2)*not(nf1-nf2)*not(nas1-nas2);
end;end;end;end;end;end;end;end;