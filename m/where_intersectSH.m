function [X,d] = where_intersectSH(r1, r2, rSH, rstart);
[X,d] = fminsearch(@(r)bragg_deviation(r,r1,r2,rSH), rstart);



function d = bragg_deviation(r, r1, r2, rSH);
k1  = (r-r1); k1 = reshape(k1,2,1); k1=k1./sqrt(k1'*k1);
k2  = (r-r2); k2 = reshape(k2,2,1); k2=k2./sqrt(k2'*k2);
kSH = (r-rSH); kSH = reshape(kSH,2,1); kSH = 2*kSH./sqrt(kSH'*kSH);

g = (k1-k2);
%d = min((2*(kSH'*g) - g'*g)^2, (2*(kSH'*g) + g'*g)^2);
d = (2*(kSH'*g) + g'*g)^2;