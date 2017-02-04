function K = normv(M)
K = zeros(size(M)); 
for i=1:size(M,2)
    %K(:,i)=(M(:,i)-min(M(:,i)))/(max((M(:,i)))-min((M(:,i))));
     K(:,i)=M(:,i)./max(abs(M(:,i)));
end;