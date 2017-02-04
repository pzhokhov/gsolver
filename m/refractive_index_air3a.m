clear;
N = 2^12;
initial_pump_central_wavelength_nm = 400;
initial_pump_central_angular_frequency_ps_minus1 = 3e5*2*pi/initial_pump_central_wavelength_nm;
spectral_width_w0 = (1.99);
spectral_window_central_angular_frequency_ps_minus1 = initial_pump_central_angular_frequency_ps_minus1;
time_discretization_step_ps = 1./(initial_pump_central_angular_frequency_ps_minus1*spectral_width_w0/2/pi);
angular_frequency_ps_minus1 = -2*pi*((((1:N) - ((N) + 2)./2)./(time_discretization_step_ps*N))) + spectral_window_central_angular_frequency_ps_minus1;
% angular_frequency_ps_minus1 = spectral_window_central_angular_frequency_ps_minus1;
wavelength_nm = 2*pi*3e5./angular_frequency_ps_minus1;
frequency_THz = angular_frequency_ps_minus1/2/pi;
c_cm_ps_minus1 = 0.03;


xc = 370;
Temperature_K = 20 + 273.15;
pressure_Pa = 101325;
h = 0.6;

air_A1_um_minus2 = 5792105;
air_A2_um_minus2 = 167917;
air_B1_um_minus2 = 238.0185;
air_B2_um_minus2 = 57.362;

c = 0.9985;
air_C1_um_minus2 = 5792105*c;
air_C2_um_minus2 = 167917*c;
air_C3_um_minus2 = 0.06*c;
air_C4_um_minus2 = 0.08*c;
air_C5_um_minus2 = 1*c;
air_C6_um_minus2 = 0.15*c;
air_D1_um_minus2 = 238.0185;
air_D2_um_minus2 = 57.362;
air_D3_um_minus2 = (1e3/4270)^2;
air_D4_um_minus2 = (1e3/2720)^2;
air_D5_um_minus2 = (1e3/50000)^2;
air_D6_um_minus2 = (1e3/6500)^2;
air_G3_um_minus2 = 0.001;
air_G4_um_minus2 = 0.0015;
air_G6_um_minus2 = 0.0005;
refractive_index_dry_air = air_A1_um_minus2./(air_B1_um_minus2 - (1/3e2*frequency_THz).^2) + air_A2_um_minus2./(air_B2_um_minus2 - (1/3e2*frequency_THz).^2);
refractive_index_dry_CO2_air = refractive_index_dry_air*(1 + 0.534e-6*(xc-450));
refractive_index_modified_dry_air = abs(real(air_C1_um_minus2./(air_D1_um_minus2 - (1/3e2*frequency_THz).^2) ...
    + air_C2_um_minus2./(air_D2_um_minus2 - (1/3e2*frequency_THz).^2) ...
    + air_C3_um_minus2./(air_D3_um_minus2 - (1/3e2*frequency_THz).^2 + 1i*air_G3_um_minus2) ...
    + air_C4_um_minus2./(air_D4_um_minus2 - (1/3e2*frequency_THz).^2 + 1i*air_G4_um_minus2) ...
    + air_C5_um_minus2./(air_D5_um_minus2 - (1/3e2*frequency_THz).^2)) ...
    + air_C6_um_minus2./(air_D6_um_minus2 - (1/3e2*frequency_THz).^2 + 1i*air_G6_um_minus2));
refractive_index_modified_dry_CO2_air = refractive_index_modified_dry_air*(1 + 0.534e-6*(xc-450));
% L_dry_CO2_air = ((refractive_index_dry_CO2_air/1e8 + 1).^2 - 1)./((refractive_index_dry_CO2_air/1e8 + 1).^2 + 2);
% L_modified_dry_CO2_air = ((refractive_index_modified_dry_CO2_air/1e8 + 1).^2 - 1)./((refractive_index_modified_dry_CO2_air/1e8 + 1).^2 + 2);

cf = 1.022;
w0 = 295.235;
w1_um_minus2 = 2.6421;
w2_um_minus4 = -0.032380;
w3_um_minus6 = 0.004028;
refractive_index_water_vapor = cf*(w0 + w1_um_minus2*(1/3e2*frequency_THz).^2 + w2_um_minus4*(1/3e2*frequency_THz).^4 + w3_um_minus6*(1/3e2*frequency_THz).^6);
% L_water_vapor = ((refractive_index_water_vapor/1e8 + 1).^2 - 1)./((refractive_index_water_vapor/1e8 + 1).^2 + 2);


xw = humidity(Temperature_K,pressure_Pa,h);
ro_a_kg_m_minus3 = density_a(Temperature_K,pressure_Pa,xw,xc);
ro_axs_kg_m_minus3 = density_axs(xc);
ro_w_kg_m_minus3 = density_w(Temperature_K,pressure_Pa,xw,xc);
ro_ws_kg_m_minus3 = density_ws(xc);

% L = ro_a_kg_m_minus3/ro_axs_kg_m_minus3*L_dry_CO2_air + ro_w_kg_m_minus3/ro_ws_kg_m_minus3*L_water_vapor;
% L_m = ro_a_kg_m_minus3/ro_axs_kg_m_minus3*L_modified_dry_CO2_air + ro_w_kg_m_minus3/ro_ws_kg_m_minus3*L_water_vapor;
% refractive_index = 1 + (sqrt((1 + 2*L)./(1 - L)) - 1);
% refractive_index_m = 1 + (sqrt((1 + 2*L_m)./(1 - L_m)) - 1);
refractive_index = (1 + 1e-8*(ro_a_kg_m_minus3/ro_axs_kg_m_minus3*refractive_index_dry_CO2_air + ro_w_kg_m_minus3/ro_ws_kg_m_minus3*refractive_index_water_vapor));
ro_a_kg_m_minus3/ro_axs_kg_m_minus3
ro_w_kg_m_minus3/ro_ws_kg_m_minus3
refractive_index_m = (1 + 1e-8*(ro_a_kg_m_minus3/ro_axs_kg_m_minus3*refractive_index_modified_dry_CO2_air + ro_w_kg_m_minus3/ro_ws_kg_m_minus3*refractive_index_water_vapor));
refractive_index_w = interp1(wavelength_nm,refractive_index_m,3900)-1
refractive_index_3w = interp1(wavelength_nm,refractive_index_m,3900/3)-1
refractive_index_5w = interp1(wavelength_nm,refractive_index_m,3900/5)-1
2*pi/(2*pi*3e5/3900*3/0.03*(refractive_index_3w - refractive_index_w))
2*pi/(2*pi*3e5/3900*5/0.03*(refractive_index_5w - refractive_index_w))

% refractive_index = (1 + 1e-8*(ro_a_kg_m_minus3/ro_axs_kg_m_minus3*refractive_index_dry_CO2_air + ro_w_kg_m_minus3/ro_ws_kg_m_minus3*refractive_index_water_vapor));
data_Mathar = load('n_Mathar.dat','-ascii');
figure(222);plot(wavelength_nm,refractive_index-1,'b' ...
  ...  ,wavelength_nm,refractive_index_dry_CO2_air/1e8,'r' ...
  ...  ,wavelength_nm,refractive_index_modified_dry_CO2_air/1e8,'g' ...
    ,data_Mathar(1,:),data_Mathar(2,:),'r' ...
    ,wavelength_nm,refractive_index_m-1,'k' ...
  );
set(gca,'xlim',[0000 6000],'ylim',[2.66e-4 2.74e-4]);

% figure(2);plot(wavelength_nm,refractive_index-1-refractive_index_dry_CO2_air/1e8,'b');
% set(gca,'xlim',[800 5500]);

% wave_number_cm_minus1 = angular_frequency_ps_minus1.*refractive_index/c_cm_ps_minus1;
% dispersion1_ps_cm_minus1 = diff(wave_number_cm_minus1)/(2*pi/(time_discretization_step_ps*N));
% dispersion1_ps_cm_minus1(N) = 0;
% dispersion_ps2_cm_minus1 = diff(wave_number_cm_minus1,2)/(2*pi/(time_discretization_step_ps*N))^2;
% dispersion_ps2_cm_minus1(N-1) = 0;dispersion_ps2_cm_minus1(N) = 0;
% dispersion3_ps3_cm_minus1 = diff(wave_number_cm_minus1,3)/(2*pi/(time_discretization_step_ps*N))^3;
% dispersion3_ps3_cm_minus1(N-2) = 0;dispersion3_ps3_cm_minus1(N-1) = 0;dispersion3_ps3_cm_minus1(N) = 0;
% figure(2);plot(wavelength_nm,dispersion_ps2_cm_minus1*100,'b',wavelength_nm,zeros(1,N),'g');
% set(gca,'xlim',[500 3900],'ylim',[0 0.0001]);
% figure(3);plot(wavelength_nm,-2*pi*c_cm_ps_minus1*dispersion_ps2_cm_minus1./wavelength_nm./(wavelength_nm*1e-12),'b',wavelength_nm,zeros(1,N),'g');
% set(gca,'xlim',[1100 1700],'ylim',[-0.030 0]);
% figure(4);plot(wavelength_nm,dispersion3_ps3_cm_minus1*100,'b',wavelength_nm,zeros(1,N),'g');
% figure(5);plot(wavelength_nm,wave_number_cm_minus1);
% set(gca,'xlim',[1300 3900]);

% savedata = [wavelength_nm;-2*pi*c_cm_ps_minus1*dispersion_ps2_cm_minus1./wavelength_nm./(wavelength_nm*1e-12)].';
% save('dispersion_air_ps_nm_minus1_km_minus1.dat','savedata','-ascii');
% 
savedata = [wavelength_nm;refractive_index_m-1;refractive_index-1].';
save('refractive_index_minus_1_air.dat','savedata','-ascii');