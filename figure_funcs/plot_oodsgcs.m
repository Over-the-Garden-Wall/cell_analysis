C = get_constants;
oocns = C.type.oodsgc;
colmap = make_colormap(13);

figure; hold on
for n = 1:length(oocns)
    cn = oocns(n);
    
    c_d = cell_data(cn);
    soma_loc = c_d.get_midpoint(true);
     
    scatter(soma_loc(2), soma_loc(3), '*', 'markerEdgeColor', [0 0 0]);
    text(soma_loc(2), soma_loc(3), [' ' num2str(cn)]);
end

set(gca, 'YLim', [0 3.5*10^5], 'XLim',[0 3.5*10^5]);
    set(gcf, 'Position', [0 0 800 800]);
    set(gca, 'Position', [.1 .1 .8 .8]);
    


% for n = 1:length(oocns)
%     cn = oocns(n);
%     
%     c_d = cell_data(cn);
%     soma_loc = c_d.get_midpoint(true);
%     figure; 
%     set(gcf, 'Position', [0 0 1400 600]);
%     
%     subplot(1,2,1); hold on
%     
%     plot_cells(cn, 1, .01, colmap(n,:));
%     scatter(soma_loc(2), soma_loc(3), '*', 'markerEdgeColor', [0 0 0]);
%     title(num2str(cn));
%     
%     subplot(1,2,2); hold on
%     plot_cells(cn, 3, .01, colmap(n,:));
%     
% %     saveas(gcf, ['vis' num2str(cn) '.jpg']);
%     close all
% end

    
    