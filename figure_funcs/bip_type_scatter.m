type_num = 1;



types = {'t1','t2','t3a','t3b','t4'};
C = get_constants;

cmap = colormap('Lines');
cmap = cmap + randn(size(cmap))/10;
cmap(cmap>1) = 1;
cmap(cmap<0) = 0;


% for k = 1:length(types);
    cell_nums = C.type.(types{type_num});
    
    figure; hold all;
    
    for c = cell_nums;
        cell_dat = cell_data(c);
        p = cell_dat.get_surface;
        p = condense_points(p(:,2:3),10^2);
        meanp = mean(p);
        h = scatter(p(:,1),p(:,2),1);
%         ecol = get(h, 'MarkerEdgeColor');
%         set(h, 'MarkerFaceColor', ecol);
%         plot_cells(c,1, .1, cmap(mod(c,64)+1,:));
    end
    for c = cell_nums
        cell_dat = cell_data(c);
        p = cell_dat.get_surface;
        p = p(:,2:3);
        meanp = mean(p);
        text(meanp(1),meanp(2), num2str(c-60000));
    end
% end
        
    
    prep_figure(gcf,gca);
    set(gca, 'XLim', [.4*10^5 2.4*10^5])
    set(gca, 'YLim', [.4*10^5 2.4*10^5])
    