function [h, cb, lA_, x_, y_]=imagesc_(x,y,A)
% function imagesc_(A)
% Color plot with interpolation.

if (nargin == 1)
    A = x;
    x = 1:size(A,2);
    y = 1:size(A,1);
end;
if (nargin == 2)
   A = y;
   y = 1:size(A,1);
end;

lA = A; 
xislin = false; yislin = false; 
if (diff(diff(x))==0) xislin=true; end;
if (diff(diff(y))==0) yislin=true; end;

if (xislin && yislin)
    x_ = x;
    y_ = y;
    lA_ = lA;
else
    x_ = min(x):(max(x)-min(x))/(4*length(x)):max(x); x_=x_';
    y_ = min(y):(max(y)-min(y))/(4*length(y)):max(y);
    lA_ = interp2(x,y,lA, x_, y_, 'nearest'); 
end;

if (nargout > 0) 
  h = imagesc(x_, y_, lA_); 
else 
    imagesc(x_, y_, lA_); 
end;
set(gca, 'YDir', 'normal');
set(gca, 'CLim', [min(A(:)), max(A(:))]);
cb = colorbar; 


