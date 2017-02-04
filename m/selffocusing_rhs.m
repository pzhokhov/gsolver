function r = selffocusing_rhs(x,y, beta)
r(1) = y(2);
r(2) = -y(2) - abs(y(1))^2*y(1) - beta*abs(y(1))^4*y(1);

