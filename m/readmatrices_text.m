function M = readmatrices_text(name_template, ds, varargin)

Ns = [];
while (iscell(varargin{1})) 
    varargin = varargin{1};
end; 

for narg=1:length(varargin)
 Ns = [Ns length(varargin{narg})];
end;


N = prod(Ns); 

M=[]; 
for j=1:N; 
 template_args = [];   
 cj = j-1;

 for narg=1:length(varargin)
    template_args = [template_args, varargin{narg}(mod(cj,Ns(narg))+1)];
    cj = floor(cj/Ns(narg)); 
end;

 
  cname_temp = sprintf(name_template, template_args); 
  cname = cname_temp;
  
%  cname = ls(cname_temp); 
% 
%  if (isempty(cname))
%      error('File %s not found!', cname_temp); 
%  end;
%   
%  while (isspace(cname(1))) cname=cname(2:end); end; 
%  for i=1:length(cname); if (cname(i)==char(10)) cname=cname(1:i-1); break; end; end;
%  
%   while (isspace(cname(end))) 
%      if (length(cname) < 2)
%          error('Something wrong with cname');
%      end;
%      cname = cname(1:end-1); 
%  end; 
%  
 fid = fopen(cname, 'r'); 
 
 if (fid < 0) 
      error('Problem with opening file %s', cname);
 end;

 D = fscanf(fid, '%e'); 
 if (j==1) 
     d1 = ds(1); 
     if (length(ds)==1) 
         d2 = length(D)/d1;
     else
         d2 = ds(2); 
     end;
     
     M = zeros([d1, d2, Ns]); 
 end;    
 
msize = d1*d2;
if (mod(length(D), msize)~=0) error('number of elements in the file %s should be a multiply of a number of elements in the first file!', cname); end;
b = length(D)/msize; 
M(:,:,j)=reshape(D(1:b:end), d1, d2);

 
fclose(fid); 
end;
