function [d,dnew] = AnalyticalApprox3(t,newEt,GP,mu,energy,filter)
h = 6.62606957e-34;                             % Planck;
eV = 1.60217646e-19;                            % Electron volt;
n = numel(newEt);
d = zeros(numel(GP),n);
dnew = zeros(numel(GP),n);
dt = t(2) - t(1);
omega = 2*pi*energy*eV/h;
for k = 1:numel(GP)                         % Global phase calculation;
    k
    clear Et;
    Et = imag(hilbert(newEt)*exp(1i*pi*(0.5+GP(k))));
    %Et(abs(Et)<1e-3*max(Et)) = 0;
    for m = 1:numel(energy)                 % Energy levels scan;
        clear omegak omegaeff domegak;
        %omegak = mu(m)*Et/h*2*pi*8.4783e-30;
        omegak = Et;
        omegaeff = sqrt(omega(m).^2 + 4*omegak.^2);
        domegak = diff(omegak)./dt*1e15;
        domegak = [domegak domegak(end)];
        f3 = omegaeff*dt*1e-15;
        f1 = zeros(1,j);
        for l = 1:j                 % t';
                f2 = zeros(1,l);
                for g = 1:l         % t'';
                    f2(g) = sum(f3(g:l));
                end
                f1(l) = sum(domegak(1:l).*sin(f2)*dt*1e-15);
        end
        d(k,:) = cumsum(2*mu(m)./omegaeff.*f1*dt*1e-15);  %int t'
    end
   dnew(k,:) = d(k,:);%.*filter;
end