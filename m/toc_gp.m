function toc_gp(p)

global graph_progress_handle; 
global graph_progress_tic;



if (isempty(graph_progress_handle))
    return;
%  if (p == 1)
%      return;
%  else
%      error('graph_progress_handle not initialized! Use tic_graphprogress first!'); 
%  end;
end;
et = toc(graph_progress_tic); 


waitbar(p, graph_progress_handle, ...
 sprintf('%g %%, elapsed %g s, left %d s, total %d s', p*100, et, round(et*(1/p-1)), round(et/p)));

if (p==1)
    pause(0.5);
    close(graph_progress_handle); 
    graph_progress_handle = []; 
end;