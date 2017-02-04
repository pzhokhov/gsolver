function xw = humidity(T_K,p_Pa,h)

% A = 1.2378847e-5;
% B = -1.9121316e-2;
% C = 33.93711047;
% D = -6.3431645e3;
% svp_Pa = exp(A*T_K^2 + B*T_K + C + D/T_K);

K1 = 1.16705214528E+03;
K2 = -7.24213167032E+05;
K3 = -1.70738469401E+01;
K4 = 1.20208247025E+04;
K5 = -3.23255503223E+06;
K6 = 1.49151086135E+01;
K7 = -4.82326573616E+03;
K8 = 4.05113405421E+05;
K9 = -2.38555575678E-01;
K10 = 6.50175348448E+02;
omega = T_K + K9/(T_K - K10);
A = omega^2 + K1*omega + K2;
B = K3*omega^2 + K4*omega + K5;
C = K6*omega^2 + K7*omega + K8;
X = -B + sqrt(B^2 - 4*A*C);
svp_Pa = 1e6*(2*C/X)^4;

alfa = 1.00062;
beta = 3.14e-8;
gamma = 5.6e-7;
t_C = T_K - 273.15;
f = alfa + beta*p_Pa + gamma*t_C^2;
xw = f*h*svp_Pa/p_Pa;
pressure = h*svp_Pa