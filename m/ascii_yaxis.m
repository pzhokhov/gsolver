function scr_plane = ascii_yaxis(plane, offsety, miny, maxy, N)

if (nargin < 4)
    error('Too few arguements!');
end;
if (nargin < 5)
    N = 5;
end;

S = size(plane,1); 
 
template = '%+3.1e --';
for i = 1:N;
    scr_y = floor((i-1)/N*(S-sum(offsety)))+offsety(1)+1; 
        y = (i-1)/N*(maxy-miny) + miny;
    str = sprintf(template, y);
    plane(S-scr_y+1,1:length(str)) = str; 
end;    

scr_plane = plane;
