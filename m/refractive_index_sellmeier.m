function n = refractive_index_sellmeier(filename, lambda)

fid = fopen(filename, 'rb');
ftype = fread(fid, 1, 'int');
N = fread(fid, 1, 'int'); 

li = fread(fid, N, 'double');
Ai = fread(fid, N, 'double');

fclose(fid);

n2 = ones(size(lambda));
for i=1:N;
    n2 = n2 + Ai(i)./(1 - (li(i)./lambda).^2); 
end;

n = sqrt(n2);