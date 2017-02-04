function [F, Nt, tmin, tmax, xnet, ynet, znet, Vg] = load_seriessectiony(mask,I,ny)



if (nargin < 2)
    I = 0;
end;    

name = deblank(ls(sprintf(mask, I(1)))); 
[F0, Nt, tmin, tmax, xnet, ynet, znet, Vg] = load_startsectiony_old(name, 1);

if (nargin < 3) 
    ny = length(ynet)/2;
end;

if (nargin < 2)
    I = 1 : length(znet);
end;

F = zeros(Nt, length(xnet), length(I));
for i = 1:length(I);
   try	
    name = deblank(ls(sprintf(mask, I(i))));
    buf = load_startsectiony(name, ny);
   catch
 	    return;
   end;
   F(:,:,i) = buf;
end;
