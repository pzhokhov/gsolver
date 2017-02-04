function K = normh_01(M)
K = zeros(size(M)); 
for i=1:size(M,1)
    minMi = min(M(i,:));
    maxMi = max(M(i,:));
    K(i,:)=(M(i,:)-minMi)/(maxMi - minMi);
end;