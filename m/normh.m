function K = normh(M)
K = zeros(size(M)); 
for i=1:size(M,1)
    K(i,:)=M(i,:)/max(abs(M(i,:)));
end;