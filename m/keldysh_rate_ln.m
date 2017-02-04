function [W, lnW, gamma, E_amplitude] = keldysh_rate_ln(lambda, refrind, IONIZATION_POTENTIAL, lnI, emass)

LIGHT_VELOCITY		=	299792458;
VACUUM_PERMITTIVITY	=	8.854187817e-12;
VACUUM_PERMEABILITY	=	1.2566370615e-6;

ELECTRON_CHARGE     =	1.60217733e-19;
ELECTRON_MASS		=	9.1093897e-31;
PLANCK_CONSTANT		=	6.6260755e-34;
PLANCK_CONSTANT_REDUCED  =	PLANCK_CONSTANT/2/pi;
VACUUM_IMPEDANCE = sqrt(VACUUM_PERMEABILITY/VACUUM_PERMITTIVITY);
    
    OMEGA0  = 2*pi*LIGHT_VELOCITY/lambda; 
    
    
	E_amplitude = exp(lnI/2)*sqrt(2/refrind*VACUUM_IMPEDANCE);
    
    if (nargin < 5) emass = 1; end; 
    
	m = emass*ELECTRON_MASS;
	gamma = OMEGA0*sqrt(m*IONIZATION_POTENTIAL)/ELECTRON_CHARGE./E_amplitude;
	Gamma = gamma.*gamma./(1+gamma.*gamma);
	Xi    = 1./(1+gamma.*gamma);
	[Kg,Eg] = ellipke(Gamma); 
	[Kx,Ex] = ellipke(Xi);

    alph = pi.*(Kg-Eg)./Ex; bet = pi*pi./4./Kx./Ex;
	
	x     = 2/pi*IONIZATION_POTENTIAL/PLANCK_CONSTANT_REDUCED/OMEGA0.*Ex./sqrt(Gamma); 
	nu    = floor(x+1) - x;
	sumN = ceil(15./alph);
    
	Q = 0;
    if (gamma > 0.05)
     for n=0:sumN
         Q = Q+exp(-n.*alph).*dawson(sqrt(bet.*(n+2.*nu)));
     end;
    else
     Q = sqrt(3*pi)./4./gamma.^2;   
    end; 
    
%    if (g>0.1)
%     sumn = ceil(10/alpha);
%     if (progress) hw=waitbar(0,'calculating...'); end;
%     for k=nu0:nu0+sumn;
%         if (progress) 
%             waitbar((k-nu0)/sumn,hw); 
%         end;
%         q = q + exp(-alpha.*(k-nu)).*dawson(sqrt(beta.*(k-nu)));
%     end;
%     else
%        q = sqrt(3*pi)./4./g.^2;
%     end; 


    Q = Q.*sqrt(pi/2./Kx);
    
	W = 2.*OMEGA0/9/pi.*(OMEGA0.*m./PLANCK_CONSTANT_REDUCED./sqrt(Gamma)).^1.5.*Q.*exp(-alph.*floor(x+1));
    if (gamma > 1e4) W = 0; end;
    
	lnW = log(2.*OMEGA0/9/pi) + 1.5.*log(OMEGA0.*m./PLANCK_CONSTANT_REDUCED./sqrt(Gamma))+log(Q)-alph.*floor(x+1);

