h=waitbar(0, 'Calculating, be patient...'); 

d12 = 0.2:0.2:20; 
Nlev = 7;
Nd = length(d12); 
di_all = zeros(4096, 25, Nlev, Nd);
tic;
for nlev=1:Nlev; 
for nd=1:Nd; 
  
  di_all(:,:,nlev, nd)=atto_qs_(cep*pi, 1e18, nulev(nlev), d12(nd),  Qs(nlev)); 
  et = toc; 
  pbar = ((nd-1)/Nlev/Nd + (nlev-1)/length(nulev));
  waitbar(pbar, h, sprintf('Time left: %g s', et*(Nd*Nlev*(1-pbar))/(nd + Nlev*(nlev-1)) )); 
end;
end; 

fdi_all = fft(di_all); 

save di_all di_all fdi_all


close(h); 