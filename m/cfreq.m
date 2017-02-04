function w = cfreq(tspan, Nt)

if (nargin == 1 && length(tspan)>1)
    Nt = length(tspan);
    tspan = (max(tspan)-min(tspan));
end;

w = -pi*Nt/tspan : 2*pi/tspan : (pi*Nt/tspan - 2*pi/tspan);

