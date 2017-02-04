function D = load_UDFstep(filename, nz)

fname = ls(filename); fname = strtrim(fname);
fid = fopen(fname,'rb');
dumpstandard_name = 'UDF1.0';
dumpstandard_name_len = 6;
dumpstandard = fread(fid, 7, 'char=>char');
if (~strncmp(dumpstandard', dumpstandard_name,dumpstandard_name_len)) error('%s - unsupported dump standard!',dumpstandard); end;
float_size = fread(fid, 1, 'int');
if (float_size ~= 8 && float_size ~= 4) error('%d - unknown float size!', float_size); end;
iscomplex  = fread(fid, 1, 'int');
if (iscomplex ~= 1 && iscomplex ~= 0) error('%d - unknown iscomplex value (must be 0 or 1)!', float_size); end;
tN         = fread(fid, 1, 'int');
xN         = fread(fid, 1, 'int');
yN         = fread(fid, 1, 'int');
zN         = fread(fid, 1, 'int');

D = zeros(tN, xN, yN, length(nz));
if iscomplex==1
 D = complex(D);
end;

char_size = 1; intsize = 4;
header_ofs = char_size*(dumpstandard_name_len+1) + 6*intsize;

for cz = 1:length(nz)
  ofs = (1+iscomplex)*tN*xN*yN*float_size*nz(cz) + header_ofs;
  res = fseek(fid, ofs, 'bof'); if (res == -1) error('nz is too big'); end;
  if (float_size==8) buf = fread(fid, tN*xN*yN*(1+iscomplex), 'double');
  else               buf = fread(fid, tN*xN*yN*(1+iscomplex), 'float'); end;
   
  if (iscomplex==1) D(:,:,:,cz) = reshape(complex(buf(1:2:end), buf(2:2:end)),tN, xN, yN);
  else              D(:,:,:,cz) = reshape(buf,                                tN, xN, yN); end;
  
end;
D = squeeze(D);

fclose(fid);
  
