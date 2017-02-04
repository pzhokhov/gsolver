function y = hffilter(w, dw)
n = 4;
beta = 10;
wmax = max(w(:))-dw;
w_ = (w-wmax)/dw;

y = ones(size(w));
y(w_>0) = exp(-beta*(w_(w_>0).^n));

