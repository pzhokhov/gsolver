function n = ZnO_n_lambda(x)
x = x/1e-6;


n=sqrt(2.81418+0.87968*x.^2./(x.^2-0.3042.^2)-0.00711*x.^2);