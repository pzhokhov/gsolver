function ro = calculate_plasmadensity(field, tmin, tmax)

    ro = zeros(size(field(:,:)));
	N_T = size(ro, 1);
    Np  = size(ro, 2);
    
    tstep = (tmax-tmin)/N_T;
    
   
    parfor np = 1:Np
     ro_ = zeros(N_T,1);
     field_ = field(:,np);
     for nt = 2 : N_T
	
	    % Plasma density calculation in performed via Runge-Khutta method of fourth order,
	    % half-step value of A is calculated using linear interpolation
	    
	    k1 = tstep*plasma_source_function(field_(nt-1), ro_(nt-1));
	    %k2 = tstep*plasma_source_function(0.5*(field(nt-1,:)+field(nt,:)), ro(nt-1,:)+k1/2);
	    %k3 = tstep*plasma_source_function(0.5*(field(nt-1,:)+field(nt,:)), ro(nt-1,:)+k2/2);
	    k2 = tstep*plasma_source_function(field_(nt),   ro_(nt-1)+k1);

	    ro_(nt)  = ro_(nt-1) + (k1 + k2)/2;
       
     end;
     ro(:,np)=ro_; 
    end;
    ro = reshape(ro, size(field)); 
end
      
    function S = plasma_source_function(E,ro)
        const_SI;
        recombination_tau = 150e-15;
        ionization_potential = 9*1.6e-19;
        ionization_potential2 = 43; 
        collision_tau = 10e-15;
        lambda0 = 800e-9;
        n0  = 1.4533;
        wavenum0 = 2*pi/lambda0*n0;
        omega0 = 2*pi*SI.c/lambda0;
        IBS_crossection = 0;% SI.e^2/SI.m*SI.mu0/wavenum0*omega0*collision_tau/(1 + collision_tau^2*omega0^2);	
        pondero_k = SI.e^2*SI.mu0/2/SI.m/wavenum0^2/SI.c;
        
        I = abs(E).^2;
        S = -ro/recombination_tau + IBS_crossection/ionization_potential./(1+pondero_k.*I).*ro.*I + keldysh_rate_ln(lambda0, n0, ionization_potential, log(I)) + PPT_rate_ln(lambda0, n0, ionization_potential2, log(I), 2);
    end  
        
    
    