function  [E,jc, jpa, j, W, en] = read_ion(name_template, size, varargin)

nE   = [name_template '_E.dat']; 
njc  = [name_template '_jc.dat']; 
njpa = [name_template '_jpa.dat'];

E = squeeze(readmatrices_text(nE, size, varargin));
jc = squeeze(readmatrices_text(njc, size, varargin)); 
jpa = squeeze(readmatrices_text(njpa, size, varargin)); 

j = jc+jpa; 

end
