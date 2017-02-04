function B = besselsum4(X1, Y1, X2, Y2, X3, Y3, X4, Y4); 

    rho1 =  sqrt(X1.^2 + Y1.^2); phi1 = atan2(Y1, X1);
    rho2 =  sqrt(X2.^2 + Y2.^2); phi2 = atan2(Y2, X2);
    rho3 =  sqrt(X3.^2 + Y3.^2); phi3 = atan2(Y3, X3);
    rho4 =  sqrt(X4.^2 + Y4.^2); phi4 = atan2(Y4, X4); 
    
    B = besselj(0, rho1).*besselj(0, rho2).*besselj(0, rho3).*besselj(0, rho4);
    for k=1:4;
        T = besselj(k, rho1).*besselj(k, rho2).*besselj(k, rho3).*besselj(k, rho4).*2.*cos(k.*(phi1 - phi2 - phi3 + phi4)); 
        B = B+T;
    end;
    
end