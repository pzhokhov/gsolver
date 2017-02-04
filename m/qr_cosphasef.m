function A=qr_cosphasef(phN)

Q = qr_expphasef(phN);
A=0.5*(Q+Q');