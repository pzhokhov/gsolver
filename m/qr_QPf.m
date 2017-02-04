function [Q,P] = qr_QPf(phN)

A = qr_af(phN);
Q = 0.5*(A+A');
P = -0.5i*(A-A');