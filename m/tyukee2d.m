
function [B] = tyukee2d(A,tol)
 S = size(A);
 B = A;
 for j=2:S(2)-1;
 for i=2:S(1)-1;
    if (abs(B(i,j)-B(i-1,j))>tol && abs(B(i,j)-B(i+1,j))>tol && abs(B(i,j)-B(i,j-1))>tol && abs(B(i,j)-B(i,j+1))>tol)
        B(i,j)=0.25*(B(i-1,j)+B(i+1,j)+B(i,j-1)+B(i,j+1));
    end;
 end;
 end;
end