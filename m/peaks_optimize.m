

h = [1]; 
Nh = length(h); 
d12 = [3.3];
Np = length(d12);
%peaks = ones(3, Nh.^2);
%wb = waitbar(0, 'Please wait...'); 
% for j=1:Np; 
% for i=1:Nh; 
%     np = i+Nh*(j-1);
%     peaks(:,np)=[1; h(i); h(j)];
% end;end;    

tic;
clear d 
for j=1:Np; 
    d(:,:,j) = atto_bloch([ones(length(h),1), h(:)].', 1e18, 10.15, d12(j), 300); 
    %d_(:,:,j) = atto_qs_([ones(length(h),1), h(:)].', 1e18, 10.15, d12(j), 300, 1); 
    %E(:,:,j)  = atto_E_peaks([ones(length(h),1), h(:)].', I(j));
    et = toc; 
    disp(sprintf('Wait... %g percent completed, elapsed time %g s, about %g s left.', j/Np*100, et, et/j*(Np-j))); 
end;
%save datadd.mat d 

%close(wb);
