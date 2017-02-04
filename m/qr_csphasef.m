function [C,S]=qr_csphasef(phN)

Q = qr_expphasef(phN);
C=0.5*(Q+Q');
S=-0.5i*(Q-Q');