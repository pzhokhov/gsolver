function [ro,Z] = calculate_plasmadensity_gas(field, tmin, tmax, Us, P, lambda0)

    const_SI;
    
    ron = P * SI.Loschmidt_number;
    
    
    if (nargin == 3) 
        P=1;
    end;
    ro_ = zeros(size(field(:,:)));
	N_T = size(ro_, 1);
    Np  = size(ro_, 2);
    Nl = length(Us);
    
    tstep = (tmax-tmin)/N_T;
    
   % h = waitbar(0, 'Calculating, please wait...');
    for np = 1:Np
        
     cro = zeros(Nl,1);
     field_ = field(:,np);
     W0  = zeros(Nl,1);  W1=W0; 
   
     for nt = 2 : N_T
		   
       lnI1 = log(abs(field_(nt).^2)); 
      if (lnI1 > 30 )
        W0 = W1; 
        for nl=1:Nl; 
            W1(nl) = PPT_rate_ln(lambda0, 1, Us(nl), lnI1, nl); 
        end;
        k1 = tstep*W0(1)*(ron - cro(1)); k2 = tstep*W1(1)*(ron -cro(1)-k1); 
        if (k1 < 0) k1=0; end; if (k2<0) k2=0; end;
        cro(1) = cro(1) + 0.5*(k1+k2); 
        if (cro(1) > ron) cro(1)=ron; end;
        ro_(nt, np) = cro(1); 
        for nl=2:nl; 
            k1 = tstep*W0(nl)*(cro(nl-1) - cro(nl)); k2 = tstep*W1(nl)*(cro(nl-1) -cro(nl)-k1); 
            if (k1 < 0) k1=0; end; if (k2<0) k2=0; end;
            cro(nl) = cro(nl)+0.5*(k1+k2);
            if (cro(nl) > cro(nl-1)) cro(nl) = cro(nl-1); end;
            ro_(nt, np) = ro_(nt, np) + cro(nl); 
        end;     
      else
       ro_(nt, np) = ro_(nt-1, np);
      end; 
     % waitbar((nt-1 + N_T*(np-1))/Np/N_T, h);
     end;
    end;
   % close(h);
    ro = reshape(ro_, size(field)); 
    Z = ro/ron;
end
      