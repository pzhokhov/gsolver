function toc_graphprogress(p)

global graph_progress_handle; 
if (isempty(graph_progress_handle)) error('graph_progress_handle not initialized! Use tic_graphprogress first!'); 

et = toc; 


waibar(p, graph_progress_handle, ...
 printf('%g percent complete, elapsed time %g s, time left %g s, total time estimate %g s', p*100, et, et*(1/p-1), et/p)));

if (p==1)
    pause(0.5);
    close(graph_progress_handle); 
    graph_progress_handle = []; 
end;