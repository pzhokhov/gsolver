%// Constants for LME
%// 18 Sept 2004
%// Copyright 2002-2005, Calerga

%// units: SI


global SI; 
SI.faraday_constant = 9.6485309e4;
SI.gravitational_constant = 6.672659e-11;
SI.hubble_constant = 3.2e-18;
SI.ice_point = 273.15;
SI.mu0 = 1.256e-6;
SI.epsilon0 = 8.854187817e-12;

SI.molar_gaz_constant = 8.314510;
SI.molar_volume = 22.41410e-3;
SI.Loschmidt_number = 2.687e25; 

SI.e = 1.60217733e-19;
SI.m = 9.1093897e-31;
SI.mu_mass = 1.8835327e-28;
SI.n_mass = 1.6749286e-27;
SI.p_mass = 1.6726231e-27;

SI.h = 6.6260755e-34;
SI.h_ = SI.h / (2 * pi);
SI.plank_mass = 2.17671e-8;

SI.c = 299792458;
SI.sound_v = 340.29205;

SI.solar_radius = 6.9599e8;
SI.stefan_boltzmann_constant = 5.67051e-8;

SI.eV = SI.e;  %electron-volt
SI.keV = SI.eV*1e3;
SI.MeV = SI.eV*1e6;
SI.GeV = SI.eV*1e9;
SI.avogadro_number = 6.02214179e23;
SI.Ry  = 13.60569253*SI.e;

SI.bohr_radius = 5.291772085936e-11; 
SI.atomic_field = 5.142208548937811e+11;
SI.hartree_energy = 4.359748124303985e-18;

SI.atomic_dipole = SI.bohr_radius*SI.e;
SI.atomic_energy = SI.hartree_energy; 
SI.atomic_frequency = SI.atomic_energy/SI.h_; 
SI.atomic_time   = 1/SI.atomic_frequency;
SI.amu = 1.66053892e-27; 