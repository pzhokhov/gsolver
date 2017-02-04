function [E, w, phix] = xangle_interpolation(A, omega, wavenum, xmin, xmax)
%function [E, phix] = xangle_intepolation(A,omega,wavenum,xmin,xmax)
%
%This function takes field amplitude A(t,x), performs Fourier transform to
%A(w,kx) and interpolates it to eqidistant angle net A(w,phi)
%phi = arcsin(kx/k(w))
%
% omega is field circular frequency. May be fftshift'ed or not. 
% wavenum - field wavenumber in the same format with omega
% xmin - minimum value of x coordinate
% xmax - maximum value of x coordinate

if (length(omega)~=size(A,1))
    error('length(omega) should be equal to size(A,1)');
end;
if (length(wavenum)~=size(A,1))
    error('length(wavenum) should be equal to size(A,1)');
end;
if (xmin >= xmax)
    error('xmin should be less than xmax');
end;

xspan = xmax-xmin;
Nt = size(A,1);
Nx = size(A,2);
if (omega(1) > omega(Nt-1))
    w = fftshift(omega);
    k = real(fftshift(wavenum));
else
    w = omega;
    k = real(wavenum);
end;

kx = (-pi*Nx/xspan+2*pi/xspan) : 2*pi/xspan : pi*Nx/xspan;
phix = -1 : 0.002 : 1;

%A_ = zeros(size(A)); A_(:, 2869:3407) = A(:,2869:3407);
fA = abs(fftshift(fft2(A)));
E = zeros(Nt, length(phix));
for nt=1:Nt;
 
    if (w(nt)) <= 0
        continue;
    end;
    sphix = kx./k(nt);
    cx = abs(sphix) < 1;
    sphix1 = sphix(cx);
 
    E(nt,:) = interp1(asin(sphix1), fA(nt,cx), phix, 'spline')...
              *(k(nt)^2).*cos(phix); 
end;



