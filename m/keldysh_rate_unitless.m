function [W, lnW, gamma] = keldysh_rate_unitless(E, al, omega)

    
    OMEGA0  = omega; 
    
    
	E_amplitude = E; 
    
       
	m = 1/al; 
	gamma = OMEGA0/E/sqrt(al);
    
	Gamma = gamma.*gamma./(1+gamma.*gamma);
	Xi    = 1./(1+gamma.*gamma);
	[Kg,Eg] = ellipke(Gamma); 
	[Kx,Ex] = ellipke(Xi);

    alph = pi.*(Kg-Eg)./Ex; bet = pi*pi./4./Kx./Ex;
	
	x     = 2/pi/OMEGA0.*Ex./sqrt(Gamma); 
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
    
	W = 2.*OMEGA0/9/pi.*(OMEGA0.*m./sqrt(Gamma)).^1.5.*Q.*exp(-alph.*floor(x+1));
    if (gamma > 1e4) W = 0; end;
    

    lnW = log(2.*OMEGA0/9/pi) + 1.5.*log(OMEGA0.*m./sqrt(Gamma))+log(Q)-alph.*floor(x+1);

