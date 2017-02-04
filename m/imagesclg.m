function [h, lA_, x_, y_]=imagesclg(x,y,A,n)
% function imagesclg(A)
% Logarithmic color plot
% plots A over n orders of magnitude

if (nargin == 1)
    A = x;
    x = 1:size(A,2);
    y = 1:size(A,1);
    n = 5;
end;
if (nargin == 2)
   A = x;
   n = y;
   x = 1:size(A,2);
   y = 1:size(A,1);
end;
if (nargin == 3)
    n = 5;
end;


lA = log10(A); 
max_lA = max(lA(:));
min_lA = max_lA - n;

lA(lA < min_lA) = min_lA;

xislin = false; yislin = false; 
if (diff(diff(x))==0) xislin=true; end;
if (diff(diff(y))==0) yislin=true; end;

if (xislin && yislin)
    x_ = x;
    y_ = y;
    lA_ = lA;
else
    x_ = min(x):(max(x)-min(x))/(2*length(x)):max(x); x_=x_';
    y_ = min(y):(max(y)-min(y))/(2*length(y)):max(y);
    lA_ = interp2(x,y,lA, x_, y_); 
end;

if (nargout > 0) 
  h = imagesc(x_, y_, lA_); 
else 
    imagesc(x_, y_, lA_); 
end;
set(gca, 'YDir', 'normal');
grid on; 

cb = colorbar;
ticks = get(cb,'yticklabel');
for i=1:size(ticks,1); 
    ticks_{i}=sprintf('%.2g',exp(log(10)*(str2double(ticks(i,:))))); 
end;

set(cb, 'yticklabel',ticks_);

    