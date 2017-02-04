function P = ion_Pfunction_3d(t, Ax, Ay, Az); 


dt = t(2)-t(1);
Nt = length(t); 
cAppp = cumsum(cos(Ax+Ay+Az))*dt/4;
sAppp = cumsum(sin(Ax+Ay+Az))*dt/4;

cAppm = cumsum(cos(Ax+Ay-Az))*dt/4;
sAppm = cumsum(sin(Ax+Ay-Az))*dt/4;

cApmp = cumsum(cos(Ax-Ay+Az))*dt/4;
sApmp = cumsum(sin(Ax-Ay+Az))*dt/4;

cApmm = cumsum(cos(Ax-Ay-Az))*dt/4;
sApmm = cumsum(sin(Ax-Ay-Az))*dt/4;


P = zeros(Nt); 

parfor nt = 1:Nt; 
    
    Cppp = cAppp - cAppp(nt);
    Sppp = sAppp - sAppp(nt);
    
    Cpmp = cApmp - cApmp(nt);
    Spmp = sApmp - sApmp(nt);
    
    Cppm = cAppm - cAppm(nt);
    Sppm = sAppm - sAppm(nt);
    
    Cpmm = cApmm - cApmm(nt);
    Spmm = sApmm - sApmm(nt);
    
    P(:,nt) = besselsum4(Cppp, Sppp, Cpmp, Spmp, Cppm, Sppm, Cpmm, Spmm); 
end; 


