c = 3*10^8;
T = 100*10^(-15); % want half of whatever 
% T window you used in the config file.
N = 512; %Put in whatever N_R you used here.
lambda = 800*10^(-9); %again change to whatever you used.

wl = zeros(1, N);
omega = 2*pi*c/lambda;

f = linspace(0, 2*omega, N);
% can leave this here if all you need is frequency.
% for wavelength units:

for i = 1:N-1;
wl(1, i+1) = 2 * pi * c/f(1, i + 1);
end

wl(1, 1) = 2e-4; %just some arbitrarily large value.