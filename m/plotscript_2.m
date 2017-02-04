nas = [10:2:20]; 

fontsize = 20;
nz_w = [2049, 2900];
for i=1:length(nas); 
    figure;
 %   subplot(length(nas)/2,2,i);
    imagesc_(t,z(sz), Ew(sz,:,nas(i)));
    ax = gca; 
    set(ax, 'YaxisLocation', 'right', 'FontSize', fontsize); 
    xlabel('time (cycles)'); ylabel('z (wavelengths)'); 
     
%      newax = axes('position', get(ax,'Position')); 
%      plot(Ww(sz, end, nas(i))/rho_c, z(sz), '-.','LineWidth', 3, 'Color', 'magenta'); 
%      set(newax, 'Color', 'none', 'XDir', 'reverse', 'Ylim', [min(z(sz)) max(z(sz))], 'XAxisLocation' ,'top', 'YAxisLocation', 'right'); 
     
    newax2 = axes('position', get(ax,'Position')); 
    l2 = plot(t, (Ww(nz_w,:,nas(i)))/rho_c, 'LineWidth', 3); 
    
    set(l2(1), 'LineStyle', '--'); 
    set(l2(2), 'LineStyle', '-.'); 
    
    
    cylim = get(newax2, 'Ylim'); ymax = cylim(2); 
    set(newax2, 'Color', 'none', 'Xlim', get(ax, 'XLim'), 'YLim', cylim,  'YAxisLocation' ,'left'); 
    set(newax2, 'FontSize', fontsize); 
    ylabel('CB population (\rho_c)'); 
    
    set(newax2, 'Position', get(ax, 'Position'));
    
%    ax_(i) = ax; 
%    newax2_(i) = newax2;
    
% 
%     
%     subplot(length(nas),2,2*i); 
%     lhs= plot(t, E(2448,:,nas(i)), t(st), E(1648,st,nas(i)), '--',  t, W(2448,:,nas(i))/rho_c, '-.', t+2*(z(2048)-z(1648)), E(1648, :, nas(i)),'LineWidth', 2); 
%     set(gca, 'Xlim', [5 12]); set(lhs(end), 'LineWidth', 1, 'Color', 'k'); 

    %set(gcf, 'PaperPositionMode','auto');
    print(gcf, sprintf('~/work/papers/breakdown/fdtd_Ew_a%g.eps', a(nas(i))), '-depsc');
end;


% for i=1:length(nas); 
%     set(newax2_(i), 'Position', get(ax_(i), 'Position'));
% end;
