function [romax, Imax, flux, znet, xmin, xmax] = parse_sectiondump(filename)

plasmaname = [filename,'.plasma'];
[A, tmin, tmax, xmin, xmax, znet] = load_sectiondump(filename,0);

Nz = length(znet);
[Nt,Nx] = size(A);
romax = zeros(Nx, Nz);
Imax  = romax;
flux = romax;
tstep = (tmax-tmin)/Nt;

for i=0:Nz-1;
	ro = load_sectiondump(plasmaname,i);
	A = abs(load_sectiondump(filename,i)).^2;

	romax(:,i+1) = squeeze(max(ro,[],1));
	Imax(:,i+1)  = squeeze(max(A, [],1));
	flux(:,i+1)  = squeeze(sum(A,1))*tstep;
end;

