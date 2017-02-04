function scr_plane = ascii_xaxis(plane, offsetx, minx, maxx, N)

if (nargin < 4)
    error('Too few arguements!');
end;
if (nargin < 5)
    N = 5;
end;

S = size(plane,2); 

template1 = '%3.1e';
template2 = '|'; 
for i = 1:N;
    scr_x = floor((i-1)/N*(S-sum(offsetx)))+offsetx(1)+1; 
        x = (i-1)/N*(maxx-minx) + minx;
    str = sprintf(template1, x);
    plane(size(plane,1)  ,scr_x:scr_x+length(str)-1) = str;
    plane(size(plane,1)-1,scr_x:scr_x+length(template2)-1) = template2;
end;    

scr_plane = plane; 
