function n = quartz_n_lambda(lambda)
const_SI;

B = [0.696163 0.4079426 0.897494];
lambda_r = [0.0684043 0.1162414 9.896161]*1e-6;

omega_r =  2*pi*light_v./lambda_r;
omega   =  2*pi*light_v./lambda;

n = ones(size(lambda));
for i = 1 : size(lambda,1)
for j = 1 : size(lambda,2)
    for m = 1:length(B)
        n(i,j) = n(i,j) + B(m)*(omega_r(m)^2)./((omega_r(m)^2) - omega(i,j)^2);
    end
end
end

n = n.^(0.5);