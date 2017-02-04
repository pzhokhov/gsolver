function eqistepsave(A, xnet, znet, filename, interpmethod)

if (nargin < 4)
    error('Too few arguements!');
end;

if (nargin < 5)
    interpmethod = 'nearest';
end;

   


Nx = length(xnet);
xstep = min(diff(xnet));
xneteq = (min(xnet) : xstep : max(xnet))';

Nz = length(znet);
zstep = min(diff(znet));
zneteq = (min(znet) : zstep : max(znet));

Aeq = interp2(znet, xnet, A, zneteq, xneteq, interpmethod);

Aeq_ = [0 zneteq; xneteq Aeq];

save(filename, 'Aeq_', '-ascii');


