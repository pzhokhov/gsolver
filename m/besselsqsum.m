function F = besselsqsum(k0, x);

Nsum = 1000; 
S = 0;
for n=k0:(k0+Nsum)
    S = S + besselj(k0, x).^2;
end;

F = S;