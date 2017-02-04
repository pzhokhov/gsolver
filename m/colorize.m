function A = colorize(S, frequency)


%hsvmap = colormap('hsv');
sizeS = size(S);

S = (S-min(S(:)))./(max(S(:))-min(S(:)));
A_R    = zeros([sizeS(2:end) 1]);
A_G    = A_R;
A_B    = A_R;
const_SI;
lambda = 2*pi.*light_v./frequency./1e-9;
sRGB   = spectrumRGB(lambda);
for nw = 1:size(S,1);
    
    A_R(:) = A_R(:) + sRGB(1,nw,1)*S(nw,:)';
    A_G(:) = A_G(:) + sRGB(1,nw,2)*S(nw,:)';
    A_B(:) = A_B(:) + sRGB(1,nw,3)*S(nw,:)';
end;

Mr = max(A_R(:)); Mg = max(A_G(:)); Mb = max(A_B(:));
mr = min(A_R(:)); mg = min(A_G(:)); mb = min(A_B(:));

M = max([Mr;Mg;Mb]);
m = min([mr;mg;mb]);

A_R = (A_R-m)./(M-m);
A_G = (A_G-m)./(M-m);
A_B = (A_B-m)./(M-m);

A = cat(length(size(S)), A_R, A_G, A_B);

    