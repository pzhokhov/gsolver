function U = group_velocity_lambda(l, f)

const_SI; 
U = zeros(size(l)); 
w = 2*pi*SI.c./l(:);
delta = 1e-4; 

w_ = [w*(1-delta) w*(1+delta)];

l_ = 2*pi*SI.c./w_;
n_ = f(l_); 

k_ = w_.*n_./SI.c; 

U(:) = diff(w_, 1, 2)./diff(k_, 1, 2); 
