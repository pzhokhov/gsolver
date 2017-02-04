function toc_dispprogress(p, st)

if (nargin >= 2)
    et = toc(st); 
else
    et = toc; 
end;
disp(sprintf('%g percent complete, elapsed time %g s, time left %g s, total time estimate %g s', p*100, et, et*(1/p-1), et/p)); 