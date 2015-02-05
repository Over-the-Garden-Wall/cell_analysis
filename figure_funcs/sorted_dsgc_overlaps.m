cell_nums = {[90001 25005 20213 20220 17080], [90002], [17161 20233], [20239 20210 20254 20245]};
lbls = {'red', 'green', 'yellow', 'purple', 'unknown'};

colors = [1 0 0; 0 .7 0; .8 .8 0; 1 0 1];
figure; hold on
for t = 1:length(cell_nums)
    
%     
%     visualize_overlaps(cell_nums{t}+1000);
%     title([lbls{t} ' OFF']);
%     saveas(gcf, ['~/data/stratification/images/overlap_' lbls{t} 'off.png']);
%     
%     visualize_overlaps(cell_nums{t}+2000);
%     title([lbls{t} ' ON']);
%     saveas(gcf, ['~/data/stratification/images/overlap_' lbls{t} 'on.png']);
%     
%     
%     plot_individual_strats(cell_nums{t}, [0 100]);
%     saveas(gcf, ['~/data/stratification/images/strats_' lbls{t} '.png']);
%     
%     
%     types = {'t1', 't2', 't3a', 't3b', 't4', 't5w', 't5h', 't5l', 'xbc', 't6', 't7', 't89', 'tRBC'};
%     cns = [];
%     for n = 1:length(types)
%         cns{n} = C.type.(types{n});
%     end
%     connectivity_scatter(cell_nums{t}, cns);
%     set(gca, 'XTick', 1:length(types), 'XTickLabel', types);
%     saveas(gcf, ['~/data/stratification/images/grossbip_' lbls{t} '.png']);
    for k = 1:length(cell_nums{t});
        c_d = cell_data(cell_nums{t}(k)+1000);
        off_mid = c_d.get_midpoint(true);
        c_d = cell_data(cell_nums{t}(k)+2000);
        on_mid = c_d.get_midpoint(true);
        
        plot([0; on_mid(2) - off_mid(2)], [0; on_mid(3) - off_mid(3)], 'lineWidth', 2, 'Color', colors(t,:));
    end
    
end


% close all
    