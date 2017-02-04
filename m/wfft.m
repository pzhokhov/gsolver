function M = wfft(x, wl); 

t = (1:length(x)).';
M = zeros(length(x));
for i=1:length(x); 
    xw = x(:).*exp(-((t-i)/wl).^2); 
    M(:,i)=fftshift(fft(xw));
end;

