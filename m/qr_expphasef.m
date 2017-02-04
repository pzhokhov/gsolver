function A=qr_expphasef(phN)

A = (qr_Nf(phN)+eye(phN^3))^0.5*qr_af(phN);