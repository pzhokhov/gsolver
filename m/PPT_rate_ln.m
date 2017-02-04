function [W, lnW, Q, g] = PPT_rate_ln(lambda, refrind, IONIZATION_POTENTIAL, lnI, Z, progress)


if (nargin < 5) 
    Z=1;
end;

if (nargin < 6) 
    progress=false;
end;


if (lnI < log(1e10))
    W=0; lnW=-inf; Q=NaN; g=NaN;
    return;
end;    


l = 1;


p = 1*1.0129e5;  %pressure, in atmospheres

T = 293;

LIGHT_VELOCITY		=	299792458;
VACUUM_PERMITTIVITY	=	8.854187817e-12;
VACUUM_PERMEABILITY	=	1.2566370615e-6;


ELECTRON_CHARGE     =	1.60217733e-19;
ELECTRON_MASS		=	9.1093897e-31;
PLANCK_CONSTANT		=	6.6260755e-34;
PLANCK_CONSTANT_REDUCED  =	PLANCK_CONSTANT/2/pi;
VACUUM_IMPEDANCE = sqrt(VACUUM_PERMEABILITY/VACUUM_PERMITTIVITY);

IONIZATION_POTENTIAL = IONIZATION_POTENTIAL*ELECTRON_CHARGE;

    OMEGA0  = 2*pi*LIGHT_VELOCITY/lambda; 
    
    
	E = exp(lnI/2)*sqrt(2/refrind*VACUUM_IMPEDANCE);
    
    Ui = IONIZATION_POTENTIAL;
	m = ELECTRON_MASS;
	g = OMEGA0*sqrt(2*m*Ui)/ELECTRON_CHARGE./E; %keldysh parameter
    Ueff = Ui*(1+1./2./g.^2);
    Uh = ELECTRON_MASS/(4*pi*VACUUM_PERMITTIVITY/ELECTRON_CHARGE*PLANCK_CONSTANT_REDUCED/ELECTRON_CHARGE)^2; %Hartree energy
    
    nu = Ueff/PLANCK_CONSTANT_REDUCED/OMEGA0;
    nu0 = floor(1+nu);
    alpha = 2*(asinh(g)-g./sqrt(1+g.^2));
    beta = 2.*g./sqrt(1+g.^2);
    E0 = sqrt(2*Ui^3*m)/ELECTRON_CHARGE/PLANCK_CONSTANT_REDUCED;
    neff = Z/sqrt(2*Ui/Uh); 
  
    Q = 0;
    if (g>0.5)
    sumN = ceil(10/alpha);
    if (progress) hw=waitbar(0,'Calculating...'); end;
    for k=nu0:nu0+sumN;
        if (progress) 
            waitbar((k-nu0)/sumN,hw); 
        end;
        Q = Q + exp(-alpha.*(k-nu)).*dawson(sqrt(beta.*(k-nu)));
    end;
    else
       Q = sqrt(3*pi)./4./g.^2;
    end; 
    
    if (progress) close(hw); end;
    Cnl2 = (2^(2*neff))/(neff*gamma(2*neff));
    
    
    W = 4.*sqrt(2)./pi.*Cnl2.*((2.*E0./E./sqrt(1+g.^2)).^(2.*neff-3/2)).*(2*l+1)...
        .*exp(-2.*nu.*(asinh(g)-g.*sqrt(1+g.^2)./(1+2.*g.^2))).*Ui.*g.^2/(1+g.^2).*Q./PLANCK_CONSTANT_REDUCED;
   
    lnW = log(W);
    
	
