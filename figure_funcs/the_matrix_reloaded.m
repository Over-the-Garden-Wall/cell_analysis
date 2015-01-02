function the_matrix

    C = get_constants;

    types = {'t1','t2','t3a','t3b','t4'};
    type_lbl = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'};
    
    
    num_types = length(types);
    num_per_type = zeros(num_types,1);
    all_nums = [];
    
    for k = 1:num_types
        num_per_type(k) = length(C.type.(types{k}));
    end
    
    num_sacs = length(C.type.sure_off_sac);
    num_cells = length(all_nums);
    
    M = cell(num_types,1);
    
    max_M = 0;
    

    sac_nums = C.type.sure_off_sac;
    sac_nums = sac_nums(randperm(length(sac_nums)));
    
    for j = 1:num_types
        bip_nums{j} = C.type.(types{j});
        bip_nums{j} = bip_nums{j}(randperm(length(bip_nums{j})));
    end
    
    
    load(C.conn_loc);
    for k = 1:num_sacs
        sub_conns = double([conns(2:end,conns(1,:)==sac_nums(k)) conns([1 3:end],conns(2,:)==sac_nums(k))]); 
        for j = 1:num_types
            if k == 1;
                M{j} = zeros(num_sacs, num_per_type(j));
            end
            for c = 1:num_per_type(j)
                cell_no = bip_nums{j}(c);
                M{j}(k,c) = sum(sub_conns(2,sub_conns(1,:)==cell_no));
            end
        end
    end
            
    for j = 1:num_types
        M{j} = M{j}*(16.5*16.5 + 16.5*25*2)/3;
    end
    
%     for k = 1:num_sacs
%         cell_dat = cell_data(sac_nums(k));
%         for j = 1:num_types
%             if k == 1;
%                 M{j} = zeros(num_sacs, num_per_type(j));
%             end
% 
%             
%             for c = 1:num_per_type(j)
%                 cell_no = bip_nums{j}(c);
%                 if cell_dat.contact_map.isKey(cell_no);
%                     M{j}(k,c) = (1+cell_dat.contact_area( ...
%                         cell_dat.contact_map(cell_no)) * ...
%                         (16.5*16.5 + 16.5*25*2)/3);
%                 end
%             end
%             max_M = max([max_M M{j}(k,:)]);
%         end
%     end
%     
%     figure; imagesc(log(M));
    
max_M = 2*10^6;
% close all

% for k = 1:5
%     figure;imagesc(M{k}, [0 max_M]); colormap hot
% end

    figure; hold all
    
    ax_hand = zeros(num_types,1);
    
    divider_width = .0125;

    set(gcf, 'Position', [0 0 240*(sum(num_per_type)/num_sacs*(1+divider_width*(num_types-1))) 240]);
    
    last_pos = .1-divider_width;
    remain_width = .8-divider_width*(num_types-1);
    
%     extra_plot = 3;%no idea why I have to do this...
    extra_plot = 3;
    for k = [1:num_types extra_plot] 
        if k == extra_plot
            extra_pos = last_pos;
        end
            
        ax_hand(k) = subplot(1,num_types,k);
        imagesc(M{k}, [0 max_M]);        
        colormap hot;
    
        set(ax_hand(k), 'TickDir', 'out');

        set(ax_hand(k), 'FontSize', 20);

        set(ax_hand(k), 'XTick', num_per_type(k)/2)
        set(ax_hand(k), 'XTickLabel', type_lbl(k));

        if k == 1
            set(ax_hand(k), 'YTick', num_sacs/2);
            set(ax_hand(k), 'YTickLabel', {'Off SACs'});
        else
            set(ax_hand(k), 'YTick', []);
        end        
        my_width = remain_width*num_per_type(k)/sum(num_per_type);
        set(ax_hand(k), 'Position', [last_pos+divider_width .1 my_width .8]);
        last_pos = last_pos + divider_width + my_width;
        
        if k == num_types
            last_pos = extra_pos;
        end
    end
        
    
%     xlim = get(gca, 'XLim');
%     ylim = get(gca, 'YLim');        
        
    
    
%     set(gca, 'LineStyle', ':')

end