function [xmin, fmin] = gpuminsearch(fun, xstart); 

Np = length(xstart); 
step_init = 0.01.*xstart; 
h = step_init;
xcur = xstart; 
grad_estimate = Np;


while sum(abs(grad_estimate).^2) > 1e-4.*Np; 

pts = zeros(Np, 2*Np); 
for i=1:length(Np); 
    pts(:,2.*i-1)   = xcur(:); pts(i,2*i-1)  =xcur(i)-h(i); 
    pts(:,2.*i)     = xcur(:); pts(i,2*i)=xcur(i)+h(i); 
end; 

%gpts = gpuArray(pts); 

[fvals] = arrayfun(fun, pts);   

%fvals = gather(gfvals);

grad_estimate = reshape(fvals(2:2:end)-fvals(1:2:end), Np, 1)./step_init;
xcur = xcur-grad_estimate*h; 
end; 
xmin = xcur; 
fmin = fvals;
return;