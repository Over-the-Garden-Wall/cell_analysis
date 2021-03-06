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

use_mean_contact_per_cell = false;
use_contacts_for_loc = true;

C = get_constants;
use_count = false;

types = {'t1', 't2', 't3a', 't3b', 't4'};
% types = {'j'};
num_types = length(types);

% group = [1 1 2 2 2];

cell_nums = cell(5,1);
for n = 1:num_types
    cell_nums{n} = C.type.(types{n});
end

conn_data = load(C.conn_loc);
fns = fieldnames(conn_data);

% sacs = C.type.off_sac;


bin_size = 15;
bins = bin_size/2:bin_size:bin_size*8.5;

num_bins = length(bins);
y = zeros(num_bins, num_types);
county = zeros(num_bins, num_types);

% small_sac = [70081       70089       70094       70096       70097 ...
%        70098       70099       70100       70101       70102 ... 
% 70104       70105       70106       70108];


for t = C.type.sure_off_sac

    ref_cell_num = t;

cell_dat = cell_data(ref_cell_num);
p = cell_dat.get_surface;
cell_mid = cell_dat.get_midpoint(true);

for d = 1:3;
    p(:,d) = p(:,d) - cell_mid(d);
end

p = p(:,2:3);

d = sqrt(sum(p.^2,2));

% d = d/1000;
% max_d = max(d);
max_d = 100000;
% disp(max_d);

all_dist = [];
all_cont = [];
for k = 1:num_types
    
    [all_dist{k}, all_cont{k}, cont_nums] = get_axial_contacts(ref_cell_num, cell_nums{k}, [15 70], use_contacts_for_loc);
%     [all_dist{k}, all_cont{k}, is_valid] = get_mean_axial_contacts(ref_cell_num, cell_nums{k}, use_count, false, false, false, [20 40]);
    
    [all_dist{k} sort_ind] = sort(all_dist{k});
    all_cont{k} = all_cont{k}(sort_ind);
    cont_nums = cont_nums(sort_ind);
    
    all_dist{k} = all_dist{k}/max_d*100;
    

    
    for n = 1:num_bins
        valid_range = all_dist{k} > bins(n)-bin_size/2 & all_dist{k} <= bins(n)+bin_size/2;
        
        if use_mean_contact_per_cell
            num_unique_cells = length(unique(cont_nums(valid_range)));
            increment = sum(all_cont{k}(valid_range))/num_unique_cells;
            if isnan(increment)
                increment = 0;
            end
            
            warning('this code is wrong');
        else
            increment = sum(all_cont{k}(valid_range)) * 291.5 / 1000 / 1000;
        end                    
            
        y(n,k) = y(n,k)+increment;
        county(n,k) = county(n,k)+sum(valid_range);
        
            
        
        
%         ste_y{k}(n) = std(all_cont{k}(valid_range)) / sqrt(length(valid_range)-1);
        
    end
    
end




end


% y = y./county;

fig_h = figure; ax_h = gca;

hold on
for k = 1:length(types)
    line_hands(k) = plot(bins, y(:,k), 'LineWidth', 2, 'Color', C.colormap(k,:));
end
% end



prep_figure(fig_h, ax_h, ...
    'xlabel', 'Distance from soma', ...
    'ylabel', 'Total Contact (um^2)', ...
    'legend', {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'});




% figure;
% plot(bins, y(:,3)./y(:,1), 'LineWidth', 2);
%  
% 
% set(gca, 'FontSize', 16);
% set(gca, 'box', 'off');
% set(gcf, 'Position', [0 0 1000 1000]);
% title('Ratio of bipolar type 3a contacts to type 1 contacts')


% saveas(gcf, 'sac_t3a_to_t1_ratio_totalc.eps');

% legend(line_hands, {'Type 1/2 Bipolar Cells', 'Type 3/4 Bipolar Cells'});

% figure;gca; hold(gca,'all');
% for k = 1:2
%     scatter(all_dist{k}, all_cont{k}, '*');
% end
% %         for n = 1:length(all_dist)
% %
% % %             text(all_dist(n), all_cont(n) + max(all_cont)/30, num2str(cell_nums(n)));
% % %             SA = get_size_stats(cell_nums(n));
% % %             all_cont(n) = all_cont(n)/SA;
% %         end
% 
% legend({'Type 1/2 to J contacts', 'Type 3/4 to J contacts'});
% 
% 
% title( 'contact area by distance from soma scatter')
% ylabel( 'contact area (in edges)')
% xlabel( 'distance along soma-distal axis (microns)');

% 
% figure;
% plot(bins, y(:,5)./y(:,1));
% 
% figure;
% plot(bins, county, 'LineWidth', 2);


%
%     figure; hold all
%
%     for n = 1:length(cell_nums);
%         plot(all_dist{n}, all_cont{n})
%     end
%     end

% end