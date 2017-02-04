function [ro, I, flux] = load_parseddump_float(filename)
fid = fopen(filename,'rb');
if (fid < 0) 
    error('Unable to open file %s',filename);
end;

Nx = fread(fid, 1, 'int');
Nz = fread(fid, 1, 'int');

ro   = zeros(Nx,Nz);
I    = zeros(Nx,Nz);
flux = zeros(Nx,Nz);

% buf = fread(fid, 3*Nx*Nz, 'double');
% 
% for nz=1:Nz;
%  ro  (:,nz) = buf(((nz-1)*3)  *Nx + 1:((nz-1)*3+1)*Nx);
%  I   (:,nz) = buf(((nz-1)*3+1)*Nx + 1:((nz-1)*3+2)*Nx);
%  flux(:,nz) = buf(((nz-1)*3+2)*Nx + 1:((nz-1)*3+3)*Nx);
% end
% 

ro   (:) = fread(fid, Nx*Nz, 'float');
I    (:) = fread(fid, Nx*Nz, 'float');
f        = fread(fid, Nx*Nz, 'float');
flux (1:length(f)) = f;

fclose(fid);
end
    