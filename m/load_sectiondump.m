function [A, tmin, tmax, xmin, xmax, znetfull, znet, Vg, omega0, wavenum0, omega, wavenum] = load_sectiondump(filename, aN)

fid = fopen(filename, 'rb');
if (fid < 0) 
	error('Unable to open file %s', filename); 
end;

omega0 = fread(fid, 1, 'double');
wavenum0 = fread(fid, 2, 'double'); wavenum0 = complex(wavenum0(1), wavenum0(2));
Vg = fread(fid, 1, 'double');
Nt = fread(fid, 1, 'int');
tmin = fread(fid, 1, 'double');
tmax = fread(fid, 1, 'double');

omega  = fread(fid, Nt, 'double');
wavenum     = fread(fid, 2*Nt, 'double');

wavenum    = complex(   wavenum(1:2:end),    wavenum(2:2:end));

Nx = fread(fid, 1, 'int');
xmin = fread(fid, 1, 'double');
xmax = fread(fid, 1, 'double');

Ny = fread(fid, 1, 'int');
ymin = fread(fid, 1, 'double');
ymax = fread(fid, 1, 'double');

Nz   = fread(fid, 1, 'int');
znetfull = fread(fid, Nz, 'double');

chunkN = 2*Nt*Nx+1;
if (nargin < 2)
    aN = 1:Nz;
end;
N = length(aN);
B = zeros(chunkN, N);
ofs = (10+3*Nt+Nz + aN*chunkN)*8 + 16;
for cz = 1:N
    fseek(fid, ofs(cz), 'bof');
	buf = fread(fid, chunkN, 'double');
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
