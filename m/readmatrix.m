function M = readmatrix(name, N)
fid = fopen(name, 'rb'); 

if (fid < 0) 
    error('File %s not found!', name);
end;

fseek(fid, 0,'eof');
s = ftell(fid)/8; 
fseek(fid, 0, 'bof'); 

T1 = fread(fid,s,'double');
fclose(fid); 

M = reshape(T1, N, length(T1)/N); 


