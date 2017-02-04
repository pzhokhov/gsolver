function ff = iron_gauss(f, N1, N2)
% function ff = iron_gauss(f, N1, N2);
% Ironing with gaussian window
% N1 - half-width of the window in first dimension
% N2 (optional) - half-width of the window in the second dimension

ff = f;
tau = (-2*N1 : 2*N1)';
filter = exp(-(tau./N1).^2);
filter = filter/sum(filter)*ones(1,size(ff,2));

for i = 2*N1+1 : size(f,1)-2*N1
    ff(i,:) = sum(f(i-2*N1:i+2*N1, :).*filter);
end;
   
if exist('N2','var')
    f = ff; 
    tau = -2*N2 : 2*N2;
    filter = exp(-(tau./N2).^2);
    filter = ones(size(ff,1),1)*filter/sum(filter);
    
    for i = 2*N2+1 : size(f,2)-2*N2
       ff(:,i) = sum(f(:,i-2*N2:i+2*N2).*filter,2);
    end;

end;
