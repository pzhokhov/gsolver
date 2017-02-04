function M = readmatrix_text(name, varargin)
fid = fopen(name, 'r'); 
if (fid < 0) error('Problem with opening file %s', name);
end;

D = fscanf(fid, '%e'); 

for n = 1:length(varargin)
    N(n)=varargin{n};
end;
N=[N, length(D)/prod(N)]; 

M = reshape(D, N); 


fclose(fid); 