function H = hankel_transform_eqidist(F);


K = 128; 

F = F(:);

N = size(F,1); 
nF = F.*((0:(N-1))');  
fnF = fft(nF); 

for l=0:(N-1);
    H(l+1)=0; 
 for k=0:(K-1);
     g = l*cos(k*2*pi/K);
     cg = round(g); 
     if cg <0 cg=cg+N; end;
     Phikl = fnF(cg+1); % fnF(cg+1) + (fnF(cg+2)-fnF(cg+1))*(g-cg);
     H(l+1) = H(l+1)+Phikl;
 end;
 H(l+1)=H(l+1)/K;
end;
