
function [I, I0] = spatial_coherence_examine(Acs, omega, kx, dk)

load const_SI light_v;
kx = kx(3500:4500);

%phi = asin(kx./k);
phii = -0.07:0.0001:0.05;

load glasses BGG23
filter = zeros(length(omega),1);
filter(560:800) = interp1(2*pi*light_v./BGG23(:,1)./1e-9, BGG23(:,2), omega(560:800));
fA = complex(zeros(size(Acs)));
A =  fA;
for i=1:size(Acs,3); 
  fA_ = fftshift(fft2(Acs(:,:,i)));
  fA(:,:,i) = fA_.*(filter*ones(1,size(fA_,2)));
  fA(:, 1:3400,i)=0;
  %fA(:, 4500:end,i)=0;
  fA(:,3700:end,i)=0;
  A(:,:,i) = ifft2(ifftshift(fA(:,:,i))); 
  %fA(:,:,i) = Acs(:,:,i); 
end;

n0 = 6400;

nmax = 500; 
A1 = complex(zeros(size(fA,1), nmax, size(fA,3)));
A2 = A1;
A12  = A2;

for i=1:size(A,3); 
for n=0:nmax-1; 
    %nd = floor(-0.2280*n);
    d = floor(dk*n)+1;
    
    
    A1(:,n+1,i)  = A(:,n0+n,i);
    A2(:,n+1,i)  = A(:,n0-n,i);
    
    %A1(nw, n0+n+, i) = fA(nw, :, i);
    %A2(nw, n0-n+1, i) = fA(nw, :, i);
    
    A12(:,n+1,i) = cat(1, zeros(d-1,1), 2*(A1(d:end,n+1,i).*conj(A2(1:(end-d+1),n+1,i))));
end; 
end; 
mA12 = sum(sum(A12,1),3);
mA1  = sum(sum(abs(A1).^2,1),3);
mA2  = sum(sum(abs(A2).^2,1),3);

I = mA12./(mA1 + mA2);
I0 = (mA1+mA2)/2;

    