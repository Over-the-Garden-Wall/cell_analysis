function bip_to_j_contarea_X_dist
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
use_count = false;

types = {'t1', 't2', 't3a', 't3b', 't4'};
num_types = length(types);

group = [1 1 2 2 2];

cell_nums = cell(6,1);
for n = 1:num_types
    cell_nums{n} = C.type.(types{n});
end

conn_data = load(C.conn_loc);
fns = fieldnames(conn_data);

% sacs = C.type.off_sac;


bin_size = 10;
bins = 10:bin_size:90;

num_bins = length(bins);
y = zeros(num_bins, num_types);
county = zeros(num_bins, num_types);

cell_dat = cell_data(10010);
p = cell_dat.get_surface;
cell_mid = cell_dat.get_midpoint(true);

for d = 1:3;
    p(:,d) = p(:,d) - cell_mid(d);
end

p = p(:,2:3);
d = p(:,1)*cell_dat.dist_axis(1) + p(:,2)*cell_dat.dist_axis(2);
% d = d/1000;
max_d = max(d);
disp(max_d);

for k = 1:num_types
    
    [all_dist{k}, all_cont{k}, is_valid] = get_mean_axial_contacts(10010, cell_nums{k}, use_count);
    
    cell_nums{k} = cell_nums{k}(is_valid);
    all_dist{k} = all_dist{k}(is_valid);
    all_cont{k} = all_cont{k}(is_valid);
    
    %
    %
    %         cell_dist = cell(length(cell_nums),1);
    %         cell_cont = cell(length(cell_nums),1);
    %
    %         all_dist = [];
    %         all_cont = [];
    %
    %         for n = 1:length(cell_nums);
    %             for k = 1:length(sacs);
    %                 cell_dist{n} = [cell_dist{n}; dist{k}(n)];
    %                 cell_cont{n} = [cell_cont{n}; cont_area{k}(n)];
    %             end
    %
    %             [cell_dist{n} sort_ind] = sort(cell_dist{n});
    %             cell_cont{n} = cell_cont{n}(sort_ind);
    %
    %             all_dist = [all_dist; cell_dist{n}];
    %             all_cont = [all_cont; cell_cont{n}];
    %
    %         end
    [all_dist{k} sort_ind] = sort(all_dist{k});
    all_cont{k} = all_cont{k}(sort_ind);
    cell_nums{k} = cell_nums{k}(sort_ind);
    %         all_cont = all_cont;
    
    
    %     plot(all_dist, all_cont);
    
    all_dist{k} = all_dist{k}/max_d*100;
    
    %         x = 1:C.dist_bin:ceil(max(all_dist)/C.dist_bin)*C.dist_bin;
    
    
    for n = 1:num_bins
        valid_range = all_dist{k} > bins(n)-bin_size/2 & all_dist{k} <= bins(n)+bin_size/2;
        
        y(n,k) = sum(all_cont{k}(valid_range));
        county(n,k) = sum(valid_range);
%         ste_y{k}(n) = std(all_cont{k}(valid_range)) / sqrt(length(valid_range)-1);
        
    end
    
    %
    %         x = 1:C.dist_bin:max(all_dist);
    %         y = zeros(size(x));
    %         ste_y = zeros(size(x));
    %
    %         for k = 1:length(x)
    %
    %             valid_dist = all_dist>=x(k) & all_dist<x(k)+C.dist_bin;
    %             if any(valid_dist)
    %                 y(k) = mean(all_cont(valid_dist));
    %                 ste_y(k) = std(all_cont(valid_dist)) / sqrt(sum(valid_dist)-1);
    %             end
    %         end
    %         x = x + C.dist_bin/2;
    %
%     ste_y{k}(isnan(ste_y{k})) = 0;
    
end

figure;
    line_hands = plot(bins, y, 'LineWidth', 2);
% end


title( 'Connectivity of Bipolar Cells with J-Cell')
ylabel( 'Observed connectivity (edges)')
xlabel( 'Distance along directionally selective axis (microns from soma)');
legend(types);

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

figure;
plot(bins, y(:,5)./y(:,1));

figure;
plot(bins, county, 'LineWidth', 2);


%
%     figure; hold all
%
%     for n = 1:length(cell_nums);
%         plot(all_dist{n}, all_cont{n})
%     end
%     end

end