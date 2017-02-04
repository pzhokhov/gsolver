function  o = least_squares_test(x,y)

   N = length(x);
		Sx=0; Sxx=0;
		Sy=0; Sxy=0;
        
	for i=1:N
			xc =  x(i);
			yc =  y(i);
			Sx = Sx+xc;
			Sxx = Sxx+xc*xc;
			Sxy = Sxy+xc*yc;
			Sy = Sy+yc;
    end;
    Sx = sum(x);
    Sxy = sum(x.*y);
    Sxx = sum(x.^2);
    Sy  = sum(y);
		a =   (Sxy*N - Sx*Sy)/(Sxx*N - Sx*Sx);
		b =   (Sy - Sx*a)/N;
        
		o(1) =a;
		o(2) =b;
return;
