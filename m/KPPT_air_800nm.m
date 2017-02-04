clear;
dawson_accuracy = 5e-3;
effective_electron_mass_kg = double(9.1e-31);
Plank_constant_J_ps = double(1.05e-22);
electron_charge_C = double(1.6e-19);
hidrogen_ionization_potential_J = double(13.6*1.6e-19);
magnetic_constant_J_cm_minus1_A_minus2 = double(4*pi*1e-9);
c_cm_ps_minus1 = double(0.03);
angular_central_frequency_ps_minus1 = 2*pi*(3e5/800);
angular_frequency_atomic_unit_ps_minus1 = double(4.1e4);
ionization_potential_J = 12.063*1.6e-19;
nu_0 = (ionization_potential_J/Plank_constant_J_ps/angular_central_frequency_ps_minus1);
K = fix(nu_0+1)
ionization_degree = double(1);
pressure_atm = 0.2;

powerI = double(-1:0.01:2.2) +12;
N = length(powerI);
intensity_TW_cm_minus2 = 10.^powerI*1e-12;
Keldysh_parameter = angular_central_frequency_ps_minus1/electron_charge_C ...
    ./sqrt(2*intensity_TW_cm_minus2*magnetic_constant_J_cm_minus1_A_minus2*c_cm_ps_minus1) ...
    *sqrt(2*effective_electron_mass_kg*ionization_potential_J)/100;
E_TV_cm = sqrt(2*intensity_TW_cm_minus2*magnetic_constant_J_cm_minus1_A_minus2*c_cm_ps_minus1);
E_0_TV_cm = 2*E_TV_cm(1)*Keldysh_parameter(1)*nu_0; % E_0_TV_cm = E_H_TV_cm*(ionization_potential_J/hidrogen_ionization_potential_J)^(3/2);
E_H_TV_cm = E_0_TV_cm*(ionization_potential_J/hidrogen_ionization_potential_J)^(-3/2);% E_H_TV_cm = double(0.00514);
% figure(100);loglog(intensity_TW_cm_minus2,Keldysh_parameter);
effective_n = ionization_degree*(ionization_potential_J/hidrogen_ionization_potential_J)^(-1/2);
l = double(0);
m = double(0);
effective_l = effective_n - 1;
C_n_l_2 = 2^(2*effective_n)/(effective_n*gamma(effective_n+effective_l+1)*gamma(effective_n-effective_l));% C_n_l_2 = 2^(2*effective_n)/(effective_n*gamma(2*effective_n));
f_l_m = ((2*l+1)*factorial(l+abs(m)))/(2^abs(m)*factorial(abs(m))*factorial(l-abs(m)));% f_l_m = double(1);

nu = ionization_potential_J/Plank_constant_J_ps/angular_central_frequency_ps_minus1*(1+1./(2*Keldysh_parameter.^2));%ok% nu = nu_0;

beta = 2*Keldysh_parameter./sqrt(1+Keldysh_parameter.^2);%ok% beta = 2;

alfa = 2*(asinh(Keldysh_parameter) - Keldysh_parameter./sqrt(1+Keldysh_parameter.^2));%ok% alfa = 2*(log(2*Keldysh_parameter)-1);

SIGMA = zeros(1,N);
for k = K:(K + 250)
    for i = 1:N
        if k>nu(i)
            SIGMA(i) = SIGMA(i) + exp(-alfa(i)*(k-nu(i)))*dawson(sqrt(beta(i)*(k-nu(i))),dawson_accuracy);
        end
    end
end

% /*multiphoton ionization
% k = K;
% SIGMA = exp(-alfa.*(k-nu)).*dawson(sqrt(beta.*(k-nu)),dawson_accuracy);
% */multiphoton ionization

A_m = 4/sqrt(3*pi)/factorial(abs(m))*Keldysh_parameter.^2./(1+Keldysh_parameter.^2).*SIGMA;% A_m = 1;

% figure(2);semilogx(intensity_TW_cm_minus2, A_m);
% set(gca,'Ylim', [0,10]);

g = 3./(2*Keldysh_parameter).*((1+1./(2*Keldysh_parameter.^2)).*asinh(Keldysh_parameter) - sqrt(1+Keldysh_parameter.^2)./(2*Keldysh_parameter));%ok% g = 3/2./Keldysh_parameter.*(log(2*Keldysh_parameter)-1/2);

ionization_rate_s_minus1 = angular_frequency_atomic_unit_ps_minus1*sqrt(6/pi)*C_n_l_2*f_l_m ...
    *(ionization_potential_J/(2*hidrogen_ionization_potential_J)) ...
    .*(2*E_0_TV_cm./(E_TV_cm.*sqrt(1+Keldysh_parameter.^2))).^(2*effective_n-3/2) ...
    .*A_m.*exp(-2*E_0_TV_cm./(3*E_TV_cm).*g)*1e12;

tunnel_ionization_rate_s_minus1 = angular_frequency_atomic_unit_ps_minus1*sqrt(6/pi)*C_n_l_2*f_l_m ...
    *(ionization_potential_J/(2*hidrogen_ionization_potential_J)) ...
    .*(2*E_0_TV_cm./(E_TV_cm)).^(2*effective_n-3/2) ...
    .*exp(-2*E_0_TV_cm*g./(3*E_TV_cm))*1e12;

multiphoton_ionization_rate_s_minus1 = angular_frequency_atomic_unit_ps_minus1/pi/sqrt(2)*C_n_l_2*f_l_m ...
    *(ionization_potential_J/(2*hidrogen_ionization_potential_J)) ...
    *4^(2*effective_n)*nu_0.^(2*effective_n-3/2) ...
    *nu_0.^(2*K)*exp(2*K-nu_0)*dawson(sqrt(2*(K-nu_0)),dawson_accuracy).*(E_TV_cm/E_0_TV_cm).^(2*K)*1e12;

% figure(3);loglog(intensity_TW_cm_minus2, A_m.*exp(-2*E_0_TV_cm./(3*E_TV_cm).*g),'b' ...
%     ,intensity_TW_cm_minus2, 4/sqrt(3*pi)*nu_0.^(2*K)*exp(2*K-nu_0)*dawson(sqrt(2*(K-nu_0)),dawson_accuracy).*(E_TV_cm/E_0_TV_cm).^(2*K),'r');

% /* sigma calculation
% ionization_potential_J = double(12.130*1.6e-19);
% K = double(8);
% angular_frequency_atomic_unit_ps_minus1*4^(2*effective_n)/pi/sqrt(2)*C_n_l_2*f_l_m ...
%     *(ionization_potential_J/(2*hidrogen_ionization_potential_J)) ...
%     *nu_0.^(2*effective_n+2*K-3/2) ...
%     *exp(2*K-nu_0) ...
%     .*dawson(sqrt(2*(K-nu_0)),dawson_accuracy)/(E_0_TV_cm^2/2/c_cm_ps_minus1/magnetic_constant_J_cm_minus1_A_minus2)^K*10^(-(K-1)*12)
% */ sigma calculation


multiphoton_ionization_cross_section_massive = polyfit(log10(intensity_TW_cm_minus2),log10(multiphoton_ionization_rate_s_minus1),1);
multiphoton_ionization_cross_section = 10^multiphoton_ionization_cross_section_massive(2)

% data = [intensity_TW_cm_minus2; ionization_rate_s_minus1; ...
%     multiphoton_ionization_rate_s_minus1; tunnel_ionization_rate_s_minus1].';
% save('ionization_rate.dat','data','-ascII')


%--------------------------------------------------------------------------
%-------------------multiphoton--------------------------------------------

% effective_electron_mass_kg = double(9.1e-31);
% Plank_constant_J_ps = double(1.05e-22);
% electron_charge_C = double(1.6e-19);
% hidrogen_ionization_potential_J = double(13.6*1.6e-19);
% magnetic_constant_J_cm_minus1_A_minus2 = double(4*pi*1e-9);
% c_cm_ps_minus1 = double(0.03);
% angular_central_frequency_ps_minus1 = 2*pi*double(375);
% angular_frequency_atomic_unit_ps_minus1 = double(4.1e4);
% ionization_potential_J = double(21.564*1.6e-19);
% nu_0 = (ionization_potential_J/Plank_constant_J_ps/angular_central_frequency_ps_minus1);
% K = fix(nu_0+1)
% ionization_degree = double(1);

powerI = double(-12:0.01:-1) +12;
N = length(powerI);
multiphoton_intensity_TW_cm_minus2 = 10.^powerI*1e-12;
Keldysh_parameter = angular_central_frequency_ps_minus1/electron_charge_C ...
    ./sqrt(2*multiphoton_intensity_TW_cm_minus2*magnetic_constant_J_cm_minus1_A_minus2*c_cm_ps_minus1) ...
    *sqrt(2*effective_electron_mass_kg*ionization_potential_J)/100;
E_TV_cm = sqrt(2*multiphoton_intensity_TW_cm_minus2*magnetic_constant_J_cm_minus1_A_minus2*c_cm_ps_minus1);
E_0_TV_cm = 2*E_TV_cm(1)*Keldysh_parameter(1)*nu_0; % E_0_TV_cm = E_H_TV_cm*(ionization_potential_J/hidrogen_ionization_potential_J)^(3/2);
E_H_TV_cm = E_0_TV_cm*(ionization_potential_J/hidrogen_ionization_potential_J)^(-3/2);% E_H_TV_cm = double(0.00514);



multiphoton_ionization_rate_s_minus1 = angular_frequency_atomic_unit_ps_minus1/pi/sqrt(2)*C_n_l_2*f_l_m ...
    *(ionization_potential_J/(2*hidrogen_ionization_potential_J)) ...
    *4^(2*effective_n)*nu_0.^(2*effective_n-3/2) ...
    *nu_0.^(2*K)*exp(2*K-nu_0)*dawson(sqrt(2*(K-nu_0)),dawson_accuracy).*(E_TV_cm/E_0_TV_cm).^(2*K)*1e12;

% figure(2);loglog(multiphoton_intensity_TW_cm_minus2, multiphoton_ionization_rate_s_minus1);
% set(gca, 'Ylim', [1e0 1e25],'XGrid','on','YGrid','on')


% data = [multiphoton_intensity_TW_cm_minus2;multiphoton_ionization_rate_s_minus1].';
% save('multiphoton_ionization_rate.dat','data','-ascII')

%--------------------------------------------------------------------------
%---------------------tunnel-----------------------------------------------

% effective_electron_mass_kg = double(9.1e-31);
% Plank_constant_J_ps = double(1.05e-22);
% electron_charge_C = double(1.6e-19);
% hidrogen_ionization_potential_J = double(13.6*1.6e-19);
% magnetic_constant_J_cm_minus1_A_minus2 = double(4*pi*1e-9);
% c_cm_ps_minus1 = double(0.03);
% angular_central_frequency_ps_minus1 = 2*pi*double(375);
% angular_frequency_atomic_unit_ps_minus1 = double(4.1e4);
% ionization_potential_J = double(21.564*1.6e-19);
% nu_0 = (ionization_potential_J/Plank_constant_J_ps/angular_central_frequency_ps_minus1);
% K = fix(nu_0+1)
% ionization_degree = double(1);

powerI = double(2.2:0.01:4) +12;
N = length(powerI);
tunnel_intensity_TW_cm_minus2 = 10.^powerI*1e-12;
Keldysh_parameter = angular_central_frequency_ps_minus1/electron_charge_C ...
    ./sqrt(2*tunnel_intensity_TW_cm_minus2*magnetic_constant_J_cm_minus1_A_minus2*c_cm_ps_minus1) ...
    *sqrt(2*effective_electron_mass_kg*ionization_potential_J)/100;
E_TV_cm = sqrt(2*tunnel_intensity_TW_cm_minus2*magnetic_constant_J_cm_minus1_A_minus2*c_cm_ps_minus1);
E_0_TV_cm = 2*E_TV_cm(1)*Keldysh_parameter(1)*nu_0; % E_0_TV_cm = E_H_TV_cm*(ionization_potential_J/hidrogen_ionization_potential_J)^(3/2);
E_H_TV_cm = E_0_TV_cm*(ionization_potential_J/hidrogen_ionization_potential_J)^(-3/2);% E_H_TV_cm = double(0.00514);



tunnel_ionization_rate_s_minus1 = angular_frequency_atomic_unit_ps_minus1*sqrt(6/pi)*C_n_l_2*f_l_m ...
    *(ionization_potential_J/(2*hidrogen_ionization_potential_J)) ...
    .*(2*E_0_TV_cm./(E_TV_cm)).^(2*effective_n-3/2) ...
    .*exp(-2*E_0_TV_cm./(3*E_TV_cm))*1e12;


% figure(3);loglog(tunnel_intensity_TW_cm_minus2, tunnel_ionization_rate_s_minus1);
% set(gca, 'Ylim', [1e0 1e25],'XGrid','on','YGrid','on')

% data = [tunnel_intensity_TW_cm_minus2; tunnel_ionization_rate_s_minus1].';
% save('tunnel_ionization_rate.dat','data','-ascII')

figure(1);loglog(multiphoton_intensity_TW_cm_minus2, multiphoton_ionization_rate_s_minus1*pressure_atm,'b' ...
    ,intensity_TW_cm_minus2(2:length(intensity_TW_cm_minus2)-1), ionization_rate_s_minus1(2:length(intensity_TW_cm_minus2)-1)*pressure_atm,'g' ...
    ,tunnel_intensity_TW_cm_minus2, tunnel_ionization_rate_s_minus1*pressure_atm,'r');
set(gca,'xlim',[1e-5 1e4], 'Ylim', [1e-40 1e20],'XGrid','on','YGrid','on')

data = [horzcat(multiphoton_intensity_TW_cm_minus2,intensity_TW_cm_minus2(2:length(intensity_TW_cm_minus2)-1),tunnel_intensity_TW_cm_minus2); ...
    horzcat(multiphoton_ionization_rate_s_minus1*pressure_atm,ionization_rate_s_minus1(2:length(intensity_TW_cm_minus2)-1)*pressure_atm,tunnel_ionization_rate_s_minus1*pressure_atm)].';
save('ionization_rate_air_800_nm_s_minus1.dat','data','-ascII')
