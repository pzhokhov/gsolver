function density_kg_m_minus3 = density_a(T_K,p_Pa,xw,xc)
Ma_kg_mol_minus1 = 1e-3*(28.9635 + 12.011e-6*(xc - 400));
R_J_mol_minus1_K_minus1 = 8.314510;
a0 = 1.58123e-6;
a1 = -2.9331e-8;
a2 = 1.1043e-10;
b0 = 5.707e-6;
b1 = -2.051e-8;
c0 = 1.9898e-4;
c1 = -2.376e-6;
d = 1.83e-11;
e = -0.765e-8;

t_C = T_K - 273.15;
Z = 1 - (p_Pa/T_K)*(a0 + a1*t_C + a2*t_C^2 + (b0 + b1*t_C)*xw + (c0 + c1*t_C)*xw^2) + (p_Pa/T_K)^2*(d + e*xw^2);
density_kg_m_minus3 = (p_Pa*Ma_kg_mol_minus1/(Z*R_J_mol_minus1_K_minus1*T_K))*(1 - xw);