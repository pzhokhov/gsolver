
function [x0, Ximin] = fitlinesCEP(cep, M, enu, I, nulev, Qs)

tic;
const_SI;
%xstart = [9.6]; 
%[x0, Ximin] = fminsearch(@(x)(Xi2(x, M, nM, I, nulev, Qs, nd, cep)), xstart, optimset('MaxFunEvals', 500, 'MaxIter', 500)); 
x = 0.2:0.2:25; 
h=waitbar(0, 'Please wait'); 
tic; 
[~,~,~,t]=atto_qs(cep(1)*pi, I, nulev(1), 1, Qs(1),1); 
t = t*SI.atomic_time/1e-15;
const_SI; 
w = cfreq(t); 
nu = w*SI.h_./SI.e*1e15; 

xi = zeros(length(x), length(nulev));

for j=1:length(nulev); 
    
nM = nulev(1)*(1-4/Qs(1))< enu & enu < nulev(1)*(1+4/Qs(1)); 
nd = nulev(1)*(1-4/Qs(1))<  nu &  nu < nulev(1)*(1+4/Qs(1)); 



for i=1:length(x); 
    tic; 
    xi(i,j)=Xi2(x(i), M, nM, I, nulev(j), Qs(j), nd, cep); 
    et=toc; 
    pbar = ((i-1)/length(nulev)/length(x) + (j-1)/length(nulev));
    waitbar(pbar, h, sprintf('ETA = %g s', et*(length(nulev)*length(x)*(1-pbar))));
end;


[Ximin(:,j), xl]=min(xi(:,j)); 
x0(j)=x(xl);
[~, cM, cfd] =Xi2(x0(j), M, nM, I, nulev(j), Qs(j), nd, cep); 
maxcM=max(cM); mincM=min(cM); maxcfd=max(cfd); mincfd=min(cfd); 
figure; plot(cep, (cM-mincM)./(maxcM-mincM), 'g--', cep, (cfd-mincfd)./(maxcfd-mincfd), 'b-', 'LineWidth', 2);  grid on; print(gcf, '-depsc', sprintf('line%d_cepfit', j)); 
figure; plot(x, xi(:,j), '-','LineWidth', 2); hold on; plot(x0(j), Ximin(j), 'r.', 'MarkerSize', 20);  grid on; print(gcf, '-depsc', sprintf('line%d_profile', j));              

end; 

close(h);  

function [f, cM, cfd] =Xi2(x, M, nM, I, nulev, Qs, nd, cep) 

fs = abs(x(1));
d = atto_qs(cep*pi, I, nulev, fs, Qs, 1);
fd = fftshift(fft(d),1); 

cM = sum(M(nM,:),1);  maxcM=max(cM); mincM=min(cM); 
cfd = sum(abs(fd(nd,:)).^2,1); maxcfd=max(cfd); mincfd=min(cfd);   
f = sum(((cM-mincM)./(maxcM-mincM) - (cfd-mincfd)./(maxcfd-mincfd)).^2);
x
f



