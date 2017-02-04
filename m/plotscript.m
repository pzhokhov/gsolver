nass = [1 6 10]; 
nas = [1 9 15]; 


st = t < 20; 
sz = -8 < z & z < 8;
xlim = [0 8]*2;
ylim = [-4 4]*2; 
fontsize = 23;

% for i=1:length(nas); 
%     subplot(2, length(nas), 2, 2*i-1);
%     imagesc_(t(st),z(sz),Es(sz,st,nass(i))); 
%     hold on; plot(t, zeros(size(t)), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k'); hold off;
%     set(gca, 'YLim', ylim, 'XLim', xlim); 
%     
%     subplot(length(nas), 2, 2*i);
%     imagesc_(t(st),z(sz),E(sz,st,nas(i))); 
%     hold on; plot(t, zeros(size(t)), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k'); hold off;
%     set(gca, 'YLim', ylim, 'XLim', xlim); 
% %  figure; 
% %  [ax, h1, h2] = plotyy(t(s), [E_(s,nas(i))*E_SI_factor, rho(s,nas(i))*Nc], t(s), deps(s,na(i))*eps_factor-26.85+1.4533^2);
% %  set(h1(1), 'LineWdith', 1, 'Color', 'k');
% %  set(h1(2), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0 0.5 0]); 
% %  set(h2(1), 'LineWidth', 2, 'LineStyle', '-', 'Color', [0 0 1]); 
% end;  

figure; 
for i=1:length(nass); 
     subplot(2,length(nass), i)
     imagesc_(t(st),z(sz),Es(sz,st,nass(i))); 
     
     hold on; plot(t, zeros(size(t)), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k'); hold off;
     set(gca, 'YLim', ylim, 'XLim', xlim, 'FontSize', fontsize); 
end;
  
nz = 2732; 

for i=1:length(nass); 
     subplot(2,length(nass), i+length(nass))
     plot(t(st), Es(nz,st, nass(i))*E_factor_VA, t(st), Ds(nz,st,nass(i))*E_factor_VA/1.4533.^2, '--', 'LineWidth', 2);
     set(gca, 'XLim', xlim, 'FontSize', fontsize);
end;

set(gcf, 'PaperPositionMode','auto');


figure;

for i=1:length(nass); 
     subplot(2,length(nas), i)
     imagesc_(t(st),z(sz),E(sz,st,nas(i))); 
     
     hold on; plot(t, zeros(size(t)), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k'); hold off;
     set(gca, 'YLim', ylim, 'XLim', xlim, 'FontSize', fontsize); 
      
end;
  
nz = 2732; 

for i=1:length(nass); 
     subplot(2,length(nass), i+length(nass))
     [ax, h1, h2] = plotyy(t(st), [E(nz,st, nas(i))*E_factor_VA;  D(1,st,nas(i))*E_factor_VA/1.4533.^2; E(nz,st, 1)*E_factor_VA*nas(i)/nas(1)], t(st), W(nz,st, nas(i))*Nc);
     set(ax, 'XLim', xlim, 'FontSize', fontsize);
     set(h1, 'LineWidth', 2); set(h1(2), 'LineStyle', '--');
     set(h2, 'LineWIdth', 2, 'LineStyle', '-.'); set(ax(2), 'YColor', 'k'); 
end;
set(gcf, 'PaperPositionMode','auto');