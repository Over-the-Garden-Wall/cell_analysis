% function bip_to_j_contarea_X_dist
% This script takes each contanct that the cells in cell_ids make onto
% off starburst amacrine cells, finds out how far they are from the soma of
% the appropriate SAC, and plots that distance by the surface area of the
% contact.
%
% the cell_ids variable can either be a vector of cell ids or a string
%
% if cell_ids is a vector of cell ids, it will plot the stratification of those cells.
%
% if cell_ids is a string, it will check if the string is 'j', 'off_sac', 't1',
% 't2', 't3a', 't3b', or 't4'. If it is none of those, there will be an error.
% Otherwise, it will plot the cells that belong to that group. The list of
% cell types can be seen in get_constants.m, and correspond to sebastian's recent designation


% close all


C = get_constants;

types = {'t1', 't2', 't3a', 't3b', 't4'};
ref_type = 'sure_off_sac';


bins = C.sac_bins;
bins_size = bins(2)-bins(1);


% types = {'j'};
num_types = length(types);

% group = [1 1 2 2 2];

cell_nums = cell(num_types,1);
for n = 1:num_types
    cell_nums{n} = C.type.(types{n});
end



num_bins = length(bins);
hull_y = zeros(num_bins, num_types);
% county = zeros(num_bins, num_types);

% small_sac = [70081       70089       70094       70096       70097 ...
%        70098       70099       70100       70101       70102 ... 
% 70104       70105       70106       70108];



num_ref = length(C.type.(ref_type));
ref_loc = zeros(num_ref, 3);bins = bin_size/2:bin_size:bin_size*8.5;

for s = 1:num_ref
    cell_dat = cell_data(C.type.(ref_type)(s));
    ref_loc(s,:) = cell_dat.get_midpoint(true);
end
    

all_dist = cell(num_types,1);
all_cont = cell(num_types,1);

for k = 1:num_types
    num_cells = length(cell_nums{k});
    all_dist{k} = zeros(num_cells,num_ref);
    all_cont{k} = ones(num_cells,num_ref);
    
    for c = 1:num_cells
        bip_dat = cell_data(cell_nums{k}(c));
        all_cont{k}(c,:) = all_cont{k}(c,:)*bip_dat.hull_area;
        bip_loc = bip_dat.get_midpoint(false);
        
        for s = 1:num_ref;
            all_dist{k}(c,s) = sqrt((bip_loc(2)-ref_loc(s,2))^2 + (bip_loc(3)-ref_loc(s,3))^2)/1000;
        end
            
    end
    all_dist{k} = all_dist{k}(:);
    all_cont{k} = all_cont{k}(:);
    
    %     [all_dist{k}, all_cont{k}, cont_nums] = get_axial_contacts(ref_cell_num, cell_nums{k}, [15 70], use_contacts_for_loc);
% %     [all_dist{k}, all_cont{k}, is_valid] = get_mean_axial_contacts(ref_cell_num, cell_nums{k}, use_count, false, false, false, [20 40]);
%     
%     [all_dist{k} sort_ind] = sort(all_dist{k});
%     all_cont{k} = all_cont{k}(sort_ind);
%     cont_nums = cont_nums(sort_ind);
%     
%     all_dist{k} = all_dist{k}/max_d*100;
    

    
    for n = 1:num_bins
        valid_range = all_dist{k} > bins(n)-bin_size/2 & all_dist{k} <= bins(n)+bin_size/2;
        
            increment = sum(all_cont{k}(valid_range));
            
        hull_y(n,k) = increment/num_ref;
%         counthull_y(n,k) = counthull_y(n,k)+sum(valid_range);
        
    end
        
        
%         ste_y{k}(n) = std(all_cont{k}(valid_range)) / sqrt(length(valid_range)-1);
        
    
end






% y = y./county;

fig_h = figure; ax_h = gca;

hold on
for k = 1:length(types)
    line_hands(k) = plot(bins, hull_y(:,k)/1000000, 'LineWidth', 2, 'Color', C.colormap(k,:));
end
% end

prep_figure(fig_h, ax_h, ...
    'xlabel', 'Distance from soma', ...
    'ylabel', 'Total Coverage (um^2)', ...
    'legend', {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'}, ...
    'YTick', [500 1000 1500 2000]);
    


% 
% 
% 
% ylabel( 'Total Coverage')
% xlabel( 'Distance from soma');
% legend({'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'});
% % legend(types);
% 
% 
% set(ax_h, 'FontSize', 27);
% set(ax_h, 'box', 'off');
% set(fig_h, 'Position', [0 0 640 640]);
% set(ax_h, 'Position', [.1 .1 .8 .8]);
% 
% set(ax_h, 'YTick', [500 1000 1500])


% end