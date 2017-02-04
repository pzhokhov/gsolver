function [d]=integral_duration(t,I, withwaitbar)
I_ = I(:,:);
d = zeros(size(I_,2),1);
t = t(:); 

if (nargin == 2); withwaitbar = 0; end;
if (withwaitbar == 1); h=waitbar(0,'Calculating, please wait...'); end;

for i=1:size(I_,2);
   normF = trapz(t, I_(:,i)); 
   mt = trapz(t, (t.*I_(:,i)))/normF;
   t2 = trapz(t, (t.^2.*I_(:,i)))/normF;
   d(i) = 2*sqrt(t2 - mt^2); 
   
   if (withwaitbar==1); waitbar(i/size(I_,2), h); end;
end; 


S = size(I); if (length(S)==2) S = [S 1]; end;
d = reshape(d, S(2:end));

if (withwaitbar==1); close(h); end;