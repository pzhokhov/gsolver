function out = peak2center_circshift(in);

in_ = in(:,:);

out_ = zeros(size(in_)); 
out = zeros(size(in)); 

Nt = size(out,1); 

[~,nmax] = max(abs(in_)); 

for i = 1:size(in_,2); 
    out_(:,i)=circshift(in_(:,i), Nt/2-nmax(i)); 
end;


out(:,:)=out_; 
