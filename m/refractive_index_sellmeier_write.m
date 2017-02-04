function refractive_index_sellmeier_write(filename, lambdai, Ai); 

fid = fopen(filename, 'wb');
fwrite(fid, 1, 'int');
fwrite(fid, length(lambda), 'int');
fwrite(fid, lambda, 'double');
fwrite(fid, Ai, 'double');
fclose(Ai);