C = get_constants;

% types = {'t1', 't3a'};
% types = {'t1', 't2', 'sure_t3a', 't3b', 'sure_t4'};
% types = {'t1', 't2', 't3a', 't3b', 't4'};
% types = {'t2', 't3a'};
% sac_nums = C.type.sure_off_sac;
dsgc_nums = C.type.gc37;

colmap = C.colormap;
colmap = colormap('Lines');

types = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4', 'BC5t', 'BC5o', 'BC5i', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC'};
type_layer = [1 1 1 1 1 2 2 2 2 2 2 2 2 2];
num_types = length(types);

use_soma = false;

num_gcs = length(dsgc_nums);
min_thresh = 10000;
% bins = C.sac_bins;

% bins = 20:20:240;
% bins = -150:30:150;
bins = [-200 200];

inner_thresh = 70;


% oodsgc_axes = zeros(num_gcs,7);
% for n = 1:num_gcs
%     [mean_axis, pref_axes] = dsgc_preferred_direction(dsgc_nums(n));
%     oodsgc_axes(n,:) = [dsgc_nums(n), mean_axis', pref_axes(:,1)', pref_axes(:,2)'];    
% end
% save('~/data/stratification/oodsgc_pref_axes.mat', 'oodsgc_axes');

load('~/data/stratification/oodsgc_pref_axes.mat');

dsgc_angles = atan2(oodsgc_axes(:,5), oodsgc_axes(:,4));
dsgc_type = 1 + (dsgc_angles > -2) + (dsgc_angles > -.5) + ...
    (dsgc_angles > 1) - 3*(dsgc_angles > 2.5);


bin_size = bins(2)-bins(1);
num_bins = length(bins);

% figure; off_ax = gca; hold all

c = colormap('Lines');

layer_depths = [20 40; 50 70];


gc_mids = {zeros(num_gcs,3), zeros(num_gcs,3)};
for s = 1:num_gcs
    s_dat = cell_data(dsgc_nums(s));
    if use_soma
        gc_mids{1}(s,:) = s_dat.get_midpoint(true);
        gc_mids{2}(s,:) = gc_mids{1}(s,:);
    else
        p = s_dat.get_surface;
        d = C.f(p(:,1));
        gc_mids{1}(s,:) = mean(p(d > layer_depths(1,1) & d < layer_depths(1,2),:));
        gc_mids{2}(s,:) = mean(p(d > layer_depths(2,1) & d < layer_depths(2,2),:));      
    end    
    
end

% for l = 1:2


mof = ones(1,100);
for k = 1:num_types-1
    mof(k) = estimate_mosaic_overlap(C.type.(types{k}));
end
% mof(9) = 1; 

figure; main_ax = gca; hold all;
    
    cell_dist = cell(num_types,1);
    cell_rad_dist = cell(num_types,1);
    cell_lat_dist = cell(num_types,1);
    
    total_contact = cell(num_types,1);
    total_vox_in_hull = cell(num_types,1);
    
    
    dsgc_type_contact = zeros(num_gcs, num_types);
    dsgc_type_overlap = zeros(num_gcs, num_types);

    dsgc_dist_contact = zeros(num_gcs, num_types, 6);
    dsgc_dist_overlap = zeros(num_gcs, num_types, 6);
    
    
for k = 1:num_types;
    l = type_layer(k);
    
    cell_nums = C.type.(types{k});
    num_cells = length(cell_nums);
    cell_dist{k} = zeros(num_gcs,num_cells);    
    cell_rad_dist{k} = zeros(num_gcs,num_cells);
    cell_lat_dist{k} = zeros(num_gcs,num_cells);
    
    for s = 1:num_gcs
        s_dat = cell_data(dsgc_nums(s));                
        soma_loc = gc_mids{l}(s,:);
        
        lat_axis = [-oodsgc_axes(s, 6), oodsgc_axes(s, 5)];
        for h = 1:num_cells
            h_dat = cell_data(cell_nums(h));
            mid_loc = h_dat.get_midpoint(false);
            
            for d = 2:3
                cell_rad_dist{k}(s,h) = cell_rad_dist{k}(s,h) + ...
                    (soma_loc(d)-mid_loc(d))^2;
                  cell_dist{k}(s,h) = cell_dist{k}(s,h) + (mid_loc(d)-soma_loc(d)) * oodsgc_axes(s,2+d) / 1000;
            
                  cell_lat_dist{k}(s,h) = cell_dist{k}(s,h) + (mid_loc(d)-soma_loc(d)) * lat_axis(d-1) / 1000;
            end
            cell_rad_dist{k}(s,h) = sqrt(cell_rad_dist{k}(s,h))/1000;
        end
    end
        
    [total_contact{k}, total_vox_in_hull{k}] = get_contact_density_whulls(dsgc_nums, cell_nums, layer_depths(l,:));
    
    
    total_contact{k} = total_contact{k} * mof(k);
    total_vox_in_hull{k} = total_vox_in_hull{k};
    
    
    for n = 1:num_gcs;
        dsgc_type_contact(n, k) = sum( total_contact{k}(n,:));
        dsgc_type_overlap(n, k) = sum( total_vox_in_hull{k}(n,:));
        
        
        dsgc_dist_contact(n, k, 1) = sum( total_contact{k}(n,cell_dist{k}(n,:) < 0)); 
        dsgc_dist_overlap(n, k, 1) = sum( total_vox_in_hull{k}(n,cell_dist{k}(n,:) < 0)); 
        
        dsgc_dist_contact(n, k, 2) = sum( total_contact{k}(n,cell_dist{k}(n,:) > 0)); 
        dsgc_dist_overlap(n, k, 2) = sum( total_vox_in_hull{k}(n,cell_dist{k}(n,:) > 0)); 
        
        dsgc_dist_contact(n, k, 3) = sum( total_contact{k}(n,cell_rad_dist{k}(n,:) < inner_thresh)); 
        dsgc_dist_overlap(n, k, 3) = sum( total_vox_in_hull{k}(n,cell_rad_dist{k}(n,:) < inner_thresh));         

        dsgc_dist_contact(n, k, 4) = sum( total_contact{k}(n,cell_rad_dist{k}(n,:) > inner_thresh)); 
        dsgc_dist_overlap(n, k, 4) = sum( total_vox_in_hull{k}(n,cell_rad_dist{k}(n,:) > inner_thresh));         

        dsgc_dist_contact(n, k, 5) = sum( total_contact{k}(n,cell_lat_dist{k}(n,:) < 0)); 
        dsgc_dist_overlap(n, k, 5) = sum( total_vox_in_hull{k}(n,cell_lat_dist{k}(n,:) < 0));         

        dsgc_dist_contact(n, k, 6) = sum( total_contact{k}(n,cell_lat_dist{k}(n,:) > 0)); 
        dsgc_dist_overlap(n, k, 6) = sum( total_vox_in_hull{k}(n,cell_lat_dist{k}(n,:) > 0));         

    end

    
end

close all


type_mean = zeros(num_types,1);
type_ste = zeros(num_types,1);

for k = 1:num_types
    is_valid = dsgc_type_overlap(:,k) > 0;
    type_mean(k) = mean(dsgc_type_contact(is_valid,k)./dsgc_type_overlap(is_valid,k));
    type_ste(k) = std(dsgc_type_contact(is_valid,k)./dsgc_type_overlap(is_valid,k)) / sqrt(sum(is_valid)-1);
end
    
figure; 

error_dot_plot(type_mean(:)*100, type_ste(:)*100, types);

type_mean = zeros(num_types, 3);
type_ste = zeros(num_types, 3);
% for t = 1:4
    is_me = true(size(dsgc_type));    
    
    for n = 1:num_types
        
        for c = 1:3
            conts = squeeze(dsgc_dist_contact(is_me,n,c*2 + [-1 0]));
            overs = squeeze(dsgc_dist_overlap(is_me,n,c*2 + [-1 0]));
            is_valid = overs > 0;
            
            type_mean(n,c) = mean(conts(is_valid(:,2),2)./overs(is_valid(:,2),2)) - ...
                mean(conts(is_valid(:,1),1)./overs(is_valid(:,1),1));
            
            type_ste(n,c) = std(conts(is_valid(:,2),2)./overs(is_valid(:,2),2))/sqrt(sum(is_valid(:,2))-1) + ...
                std(conts(is_valid(:,1),1)./overs(is_valid(:,1),1))/sqrt(sum(is_valid(:,1))-1);
            
            
        end
        
    end
        
%     type_mean(:,t) = mean(dsgc_type_contact(is_me,:)./dsgc_type_overlap(is_me,:), 1);
%     type_ste(:,t) = std(dsgc_type_contact(is_me,:)./dsgc_type_overlap(is_me,:), 1) / sqrt(sum(is_me)-1);
% end

for c = 1:3
figure;
error_dot_plot(type_mean(:,c)*100, 2*type_ste(:,c)*100, types);
hold all
t = get(gca, 'XLim');
plot(t, [0 0], 'k:');
% h = barplot_werror(type_mean(:,c)*100, type_ste(:,c)*100);
% set(gca, 'XTickLabel', types, 'TickDir', 'out');


prep_figure(gcf, gca, 'ylabel', 'contact density', 'size', [1280 320]);
%     
% 
t = get(gca, 'YLim');
set(gca, 'YTick', -2:.1:2, 'YLim', max(abs(t)) * [-1 1]);
% set(gca, 'YTick', [0:.5:10]);

end






