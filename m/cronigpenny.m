n0 = 1.4533; 
na  = 1.47;

lambda0 = 800e-9;
lambda  = lambda0/2;
theta0  = 0.4;

mtheta = 0.1:0.001:1.5;
dn = -0.1; 
n = (na+dn)./na;
d = 152e-9;
a = 800e-9./2./1.4533./sin(theta0) - d;

hw = waitbar(0, 'Calculating, please wait...');
for i=1:length(mtheta);
theta = mtheta(i);
phi = acos(cos(theta)/n);
alpha = 2*pi*d/lambda*na*n/sin(phi);
beta  = 2*pi*a/lambda*na/sin(theta);

%M1 = 1/(sin(phi)+n*sin(theta))*[2*sin(phi), -sin(phi)+n*sin(theta); sin(phi)-n*sin(theta), 2*n*sin(theta)];

A = sin(phi)+n*sin(theta); B = sin(phi)-n*sin(theta);
M1 = 1/2/n/sin(theta)*[A -B;-B A];

P1 = [exp(1i*alpha/sin(phi)), 0; 0, exp(-1i*alpha/sin(phi))];
M2 = M1^(-1);
P2 = [exp(1i*beta/sin(theta)), 0; 0, exp(-1i*beta/sin(theta))];

M3 = (M1*P1*M2*P2)^3;
T  = det(M3)/M3(2,2);
R  = -M3(2,1)/M3(2,2);
mT(i) = abs(T).^2;
mR(i) = abs(R).^2;
if abs(mT(i)+mR(i)-1)>0.001 
    mT(i)=nan;
    mR(i)=nan;
end;
waitbar(i/length(mtheta), hw);
end;
close(hw);

