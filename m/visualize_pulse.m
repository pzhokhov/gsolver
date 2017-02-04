function visualize_pulse(prefix, id, step);

name  = sprintf('%s_id%d_Z*',        prefix, id);
namep = sprintf('%s_plasma_id%d_Z*', prefix, id);
files  = dir(name);
filesp = dir(namep);

if (nargin == 2)
    step = 1;
end;

[P, z, Nt, tmin, tmax, xnet, ynet, nodes_x, nodes_y] = load_dumped_piece(files(1).name);
tnet = tmin : (tmax-tmin)/(Nt-1) : tmax;
wnet = pi*Nt/(tmax-tmin) : -2*pi/(tmax-tmin) : -pi*(Nt-2)/(tmax-tmin);
process_N = length(xnet)*length(ynet)/(nodes_x(2)-nodes_x(1)-1)/(nodes_y(2)-nodes_y(1)-1);
for i = 1:process_N*step:length(files); 
    name  = [files(i).name(1:end-10) '*'];
    namep = [filesp(i).name(1:end-10) '*'];
    [A,z] = load_dump(name); Ap = load_dump(namep);
    nt = size(A,1)/2;
    nx = size(A,2)/2;
    ny = size(A,3)/2;
    subplot(2,2,1); plot(tnet, abs(A(:,nx,ny)));                  title(sprintf('Temporal evolution, Z = %d',z)); drawnow;
    subplot(2,2,2); plot(wnet, abs(fftshift(fft(A(:,nx,ny)))));   title('Pulse spectrum');
    subplot(2,2,3); plot(xnet, abs(A(nt,:,ny)));   title('Transverse intensity distribution');
    subplot(2,2,4); plot(tnet, real(Ap(:,nx,ny)));   title('On-axis electron density');
    drawnow;
end;


