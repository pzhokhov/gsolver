function [ Phi ] = cw_action(A, tau1, tau2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (nargin == 2)
    dt = tau1(2)-tau1(1); 
    Phi = cumtrapz(exp(1i*A*sin(tau1)))*dt; 
    
elseif (nargin == 3)
    n1 = floor((tau1)/2/pi); 
    n2 = floor((tau2)/2/pi); 
    tau1 = tau1 - 2*pi*n1;
    tau2 = tau2 - 2*pi*n2;

    dt = 1e-3; 

    tau = tau1 : dt : tau2;

    Phi = trapz(exp(1i*A*sin(tau)))*dt; 

    Phi = Phi + (n2-n1)*besselj(0,A)*2*pi; 
end;    
end

