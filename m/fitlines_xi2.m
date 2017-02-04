function [ximin] = fitlines_xi2(d12, cep, I, nulev, Qs, S0, s); 

[di_]=atto_qs_(cep*pi, 1e18, nulev, d12, Qs); 
fdi = fft(di_); 
for nlev=1:7; 
fdi(:,:,nlev)=fdi(:,:,nlev)./max(max(fdi(330:500,:,nlev)));
end;
[al, mS, ximin] = fitspectra_(S0(s,:), fdi(s,:,:));



