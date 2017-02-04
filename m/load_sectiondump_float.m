function [A, tmin, tmax, xmin, xmax, znetfull, znet, Vg, omega0, wavenum0, omega, wavenum] = load_sectiondump_float(filename, aN)

fid = fopen(filename, 'rb');
if (fid < 0) 
	error('Unable to open file %s', filename); 
end;

omega0 = fread(fid, 1, 'float');
wavenum0 = fread(fid, 2, 'float'); wavenum0 = complex(wavenum0(1), wavenum0(2));
Vg = fread(fid, 1, 'float');
Nt = fread(fid, 1, 'int');
tmin = fread(fid, 1, 'float');
tmax = fread(fid, 1, 'float');

omega  = fread(fid, Nt, 'float');
wavenum     = fread(fid, 2*Nt, 'float');

wavenum    = complex(   wavenum(1:2:end),    wavenum(2:2:end));

Nx = fread(fid, 1, 'int');
xmin = fread(fid, 1, 'float');
xmax = fread(fid, 1, 'float');

Ny = fread(fid, 1, 'int');
ymin = fread(fid, 1, 'float');
ymax = fread(fid, 1, 'float');

Nz   = fread(fid, 1, 'int');
znetfull = fread(fid, Nz, 'float');

chunkN = 2*Nt*Nx+1;
if (nargin < 2)
    aN = 1:Nz;
end;
N = length(aN);
B = zeros(chunkN, N);
ofs = (10+3*Nt+Nz + aN*chunkN)*4 + 16;
for cz = 1:N
    fseek(fid, ofs(cz), 'bof');
	buf = fread(fid, chunkN, 'float');
	if (length(buf) < chunkN || feof(fid))
		break;
	end;
	%znet(cz) = buf(1);
	%A(:,:,cz) = reshape(complex(buf(2:2:end), buf(3:2:end)), Nt, Nx);
    B(:, cz) =  buf;
end;	

fclose(fid);
znet = B(1,:);
B = complex(B(2:2:end,:), B(3:2:end,:));
A = reshape(B, Nt, Nx, N); 
