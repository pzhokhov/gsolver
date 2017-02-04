function [varargout] = asciicontour(M)



default_marker = '*+#%@$&';
default_marker_c = 1;

argpair_c = 1;
windowsize = terminal_size() - 2;



scr_plane = char(' '*ones(windowsize)); 
offsetx = [10, 2];
offsety = [4,  2];
axissize   = [windowsize(2)-sum(offsetx), windowsize(1)-sum(offsety)];


maxM = max(M(:));
minM = min(M(:));
levels = minM : (maxM-minM)/9 : maxM;

%argc = 1; 
%while (argc <= nargin)
%if (nargin == argc)
%    y{argpair_c}=varargin{argc};
%    x{argpair_c}=(1:length(y{argpair_c}))';
%    marker(argpair_c) = default_marker(default_marker_c);
%    argc = argc + 1;
%    argpair_c = argpair_c+1;
%else if (nargin > argc && ischar(varargin{argc+1}))
%      y{argpair_c}=varargin{argc};
%      x{argpair_c}=(1:length(y{argpair_c}))';
%      marker(argpair_c) = varargin{argc+1};
%      argc = argc+2;
%      argpair_c = argpair_c+1;
%    else
%      x{argpair_c}=varargin{argc};
%      y{argpair_c}=varargin{argc+1};
%      if (nargin > argc+1 && ischar(varargin{argc+2}))
%          marker(argpair_c) = varargin{argc+2};
%          argc = argc+3;
%      else
%          marker(argpair_c) = default_marker(default_marker_c); default_marker_c = default_marker_c + 1;
%          argc = argc+2;
%      end;
%       argpair_c = argpair_c+1;
%    end;
%end;
%end; 
%if (length(x)-length(y)) 
%    error('Input arrays x and y must have the same lengths'); 
%end;
%
%minx = min(x{1}); maxx = max(x{1});
%miny = min(y{1}); maxy = max(y{1});
%for j = 2:argpair_c-1; 
%    cx = x{j};
%    cy = y{j};
%    minx = min([minx; cx(:)]); maxx = max([maxx; cx(:)]);
%    miny = min([miny; cy(:)]); maxy = max([maxy; cy(:)]);
%end; 
%
%for j = 1:argpair_c-1;
%    cx = x{j};
%    cy = y{j};
%    N = length(cx);

cx = 1:size(M,2);
cy = 1:size(M,1);
N = length(cx);
minx = cx(1); maxx = cx(end);
miny = cy(1); maxy = cy(end);
marker = default_marker(1);


for nx = 1:size(M,2);
for ny = 1:size(M,1)-1;
       
    scr_x = floor((cx(nx)-minx)*axissize(1)/(maxx-minx))+1+offsetx(1);
    scr_y = floor((cy(ny)-miny)*axissize(2)/(maxy-miny))+1+offsety(1);

    for l = levels;
        if (M(ny,nx) < l && l < M(ny+1,nx))
            scr_plane(windowsize(1)-scr_y+1, scr_x) = marker(1);
        end;
    end;
end;
end;

%end

scr_plane = ascii_yaxis(scr_plane, offsety, miny, maxy, 10);
scr_plane = ascii_xaxis(scr_plane, offsetx, minx, maxx, 10); 

for j = 1:windowsize(1);
    disp(scr_plane(j,:));
end;    
    
if (nargout > 0)
	varargout{1} = scr_plane;
end;
