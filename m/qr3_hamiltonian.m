function H = qr3_hamiltonian(phN, wr, wint)

w = 1;
if (nargin < 2)
    wr = w/10;
end;
if (nargin < 3)
    wint = w/100;
end;

H = zeros(phN^4);

for nf1=0:phN-1; for nf2=0:phN-1; for ns1=0:phN-1; for ns2=0:phN-1; for np1=0:phN-1; for np2=0:phN-1; for nas1=0:phN-1; for nas2=0:phN-1;
  n1 = 1 + nf1 + ns1*phN + np1*phN^2 + nas1*phN^3;
  n2 = 1 + nf2 + ns2*phN + np2*phN^2 + nas1*phN^3;
  
  H(n1,n2) = 2*pi*((wr*nf1 + w*ns1 + (w+wr)*np1 + (w+2*wr)*nas1)*not(nf1-nf2)*not(ns1-ns2)*not(np1-np2)*not(nas1-nas2) + ...
             wint * (sqrt( nf2   * ns2   *(np2 +1))*not(nf1-nf2+1)*not(ns1-ns2+1)*not(np1-np2-1)*not(nas1-nas2)        + ...
                     sqrt( nf2   * np2   *(nas2+1))*not(nf1-nf2+1)*not(ns1-ns2)  *not(np1-np2+1)*not(nas1-nas2-1))     + ...
             wint'* (sqrt((nf2+1)*(ns2+1)* np2    )*not(nf1-nf2-1)*not(ns1-ns2-1)*not(np1-np2+1)*not(nas1-nas2)        + ...
                     sqrt((nf2+1)*(np2+1)* nas2   )*not(nf1-nf2-1)*not(ns1-ns2)  *not(np1-np2-1)*not(nas1-nas2+1)));
end;end;end;end;end;end;end;end;
