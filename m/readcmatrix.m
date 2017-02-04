function M = readcmatrix(name, varargin)
fid = fopen(name, 'rb'); 

if (fid == -1) error('File %s not found!', name); end;

fseek(fid, 0,'eof');
s = ftell(fid)/16; 
fseek(fid, 0, 'bof'); 

T1 = fread(fid,2*s,'double');
T2 = complex(T1(1:2:end), T1(2:2:end)); 
fclose(fid); 

for n = 1:length(varargin)
    N(n)=varargin{n};
end;
N=[N, length(T2)/prod(N)]; 

M = reshape(T2, N); 


