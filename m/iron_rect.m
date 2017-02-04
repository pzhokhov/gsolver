function ff = iron_rect(f, N);
% function ff = iron_rect(f, N);
% Ironing with rectangular window
% N - half-width of the window
ff = f;

for i = N+1 : size(f,1)-N
    ff(i,:) = mean(f(i-N:i+N, :));
end;
   
