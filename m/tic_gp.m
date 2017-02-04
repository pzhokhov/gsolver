function tic_gp

global graph_progress_handle; 
global graph_progress_tic; 

graph_progress_tic = tic; 

if (~isempty(graph_progress_handle) && ishandle(graph_progress_handle))
    close(graph_progress_handle);
end; 

graph_progress_handle = waitbar(0, 'Wait...'); 