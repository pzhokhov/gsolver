function uA = unlin_joint(A, n1, n2)
N = size(A,1); 
if (nargin == 1)
    n1=1;
    n2=N;
end;

if (nargin == 2)
    n2 = n1;
    n1 = 1;
end;

N_ = numel(A)/N;
uA = reshape(A,N, N_); 
dA = mean(uA(n2,:)-uA(n1,:));
for j=1:N_;
 uA(:,j) = uA(:,j) - (1:N)'/(n2-n1)*dA - uA(n1,j);
end;
uA =reshape(uA, size(A));
