function [gamma,v] = fringe_visibility_fft(A)
 %filtering
 A = tyukee2d(A,1000);
 %A = iron_gauss(A,1,1); 
 [M,s1,s2] = max2d(A);
 %A = A-0.5*(A(10,s2)+A(500,s2));
 %search for minima
%  s1m1=s1; while A(s1m1-1,s2)<A(s1m1,s2); s1m1 = s1m1-1; end;
%  s1m2=s1; while A(s1m2+1,s2)<A(s1m2,s2); s1m2 = s1m2+1; end;
%  
%  m1 = A(s1m1,s2); m2 = A(s1m2,s2);
%  m = 0.5*(m1+m2);
%  gamma = (M-m)/(M+m);
  v = 1;
 %check for fringe periodicity
 if (M > 13000) 
     v = 0;
     warning('M = %d, CCD nonlinearity may affect fringe visibility measure',M); 
 end;
 %if (abs(s1m1+s1m2-2*s1) > 5) gamma=0; end;
  
  As = A(s1-64:s1+63,s2-64:s2+63);
  perim_mean = (mean(As(1,1:end))+mean(As(end,1:end))+mean(As(1:end,1))+mean(As(1:end,end)))/4;
  As = As - perim_mean;
  fAs = abs(fftshift(fft2(As)));
  fAsup =   fAs(1:50,60:70);
  fAsdown = fAs(80:end,60:70);
  fAscenter = fAs(51:79,60:70); 
   gamma = (max(fAsup(:))+max(fAsdown(:)))/max(fAs(:));
    s_up     = sum(sum(abs(fAsup)))    ;
    s_down   = sum(sum(abs(fAsdown)))  ;
    s_center = sum(sum(abs(fAscenter)));
    
    gamma = (s_up+s_down)/s_center; 
end

function [M,s1,s2] = max2d(A)
[P,Q] = max(A);
[M,s2] = max(P);
s1 = Q(s2);
end


