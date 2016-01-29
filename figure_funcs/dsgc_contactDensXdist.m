C = get_constants;

% types = {'t1', 't3a'};
% types = {'t1', 't2', 'sure_t3a', 't3b', 'sure_t4'};
% types = {'t1', 't2', 't3a', 't3b', 't4'};
% types = {'t2', 't3a'};
% sac_nums = C.type.sure_off_sac;
dsgc_nums = C.type.oodsgc;

colmap = C.colormap;
colmap = colormap('Lines');

types = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4', 'BC5t', 'BC5o', 'BC5i', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC'};
type_layer = [1 1 1 1 1 2 2 2 2 2 2 2 2 2];

use_soma = false;

num_gcs = length(dsgc_nums);
min_thresh = 10000;
% bins = C.sac_bins;

% bins = 20:20:240;
bins = -150:30:150;


% oodsgc_axes = zeros(num_gcs,7);
% for n = 1:num_gcs
%     [mean_axis, pref_axes] = dsgc_preferred_direction(dsgc_nums(n));
%     oodsgc_axes(n,:) = [dsgc_nums(n), mean_axis', pref_axes(:,1)', pref_axes(:,2)'];    
% end
% save('~/data/stratification/oodsgc_pref_axes.mat', 'oodsgc_axes');

load('~/data/stratification/oodsgc_pref_axes.mat');



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
for k = 1:length(types)-1
    mof(k) = estimate_mosaic_overlap(C.type.(types{k}));
end
% mof(9) = 1; 

figure; main_ax = gca; hold all;
    
    cell_dist = cell(length(types),1);
    cell_rad_dist = cell(length(types),1);
    
    total_contact = cell(length(types),1);
    total_vox_in_hull = cell(length(types),1);
    
for k = 1:length(types);
    l = type_layer(k);
    
    cell_nums = C.type.(types{k});
    num_cells = length(cell_nums);
    cell_dist{k} = zeros(num_gcs,num_cells);    
    cell_rad_dist{k} = zeros(num_gcs,num_cells);    
    
    for s = 1:num_gcs
        s_dat = cell_data(dsgc_nums(s));                
        soma_loc = gc_mids{l}(s,:);
        for h = 1:num_cells
            h_dat = cell_data(cell_nums(h));
            mid_loc = h_dat.get_midpoint(false);
            
            for d = 2:3
                cell_rad_dist{k}(s,h) = cell_rad_dist{k}(s,h) + ...
                    (soma_loc(d)-mid_loc(d))^2;
                  cell_dist{k}(s,h) = cell_dist{k}(s,h) + (mid_loc(d)-soma_loc(d)) * oodsgc_axes(s,3+d) / 1000;
            end
            cell_rad_dist{k}(s,h) = sqrt(cell_dist{k}(s,h))/1000;
        end
    end
        
    [total_contact{k}, total_vox_in_hull{k}] = get_contact_density_whulls(dsgc_nums, cell_nums, layer_depths(l,:));
    
    if k == 1
        all_cell_nums = cell_nums;
        all_types = ones(size(cell_nums));
        X = [total_contact{k}];
        H = [total_vox_in_hull{k}];
        D = [cell_dist{k}];
    else
        all_cell_nums = [all_cell_nums cell_nums];
        all_types = [all_types k*ones(size(cell_nums))];
        X = [X total_contact{k}];
        H = [H total_vox_in_hull{k}];
        D = [D cell_dist{k}];
    end
    
    
    
    is_valid = total_vox_in_hull{k} > min_thresh;
    
    total_contact{k} = total_contact{k}(is_valid(:))*.2915;
    total_vox_in_hull{k} = total_vox_in_hull{k}(is_valid(:));
    cell_dist{k} = cell_dist{k}(is_valid(:));
    cell_rad_dist{k} = cell_rad_dist{k}(is_valid(:));
%     
%     
    % sac external edges per surface voxel = 2.2488, though we're not using
    % that....
    
%      mof = [1.2839    1.2945    1.1136    1.0828    1.1160];
%     mof = [1.3952    1.3325    1.1761    1.1139    1.1452];
%     mof = [1.3952    1.3325    1.9112    1.2817    1.2875];
%      mof = [1.4639    1.6033    1.7456    1.3385    1.2451];
%      mof = [1.4639    1.6033    1.5456    1.2628    1.4905];
%  mof = [1.1351    1.1268    1.0796    1.0584    1.0481];
%     mof = ones(1,5);
    
    
    
    bin_num = ceil((cell_dist{k}-bins(1))/bin_size+.5);
   density = total_contact{k}./total_vox_in_hull{k}/.410*mof(k);
    
    plot_data{k} = zeros(num_bins,1);
    plot_ste{k} =  zeros(num_bins,1);
    
%     county{k} = zeros(num_bins,1);
    
    prob_nonzero =  zeros(num_bins,1);
    plot_normfact =  zeros(num_bins,1);
    
    for n = 1:num_bins
        my_dens = density(bin_num==n);
%         county{k}(n) = sum(bin_num==n);
        plot_data{k}(n) = mean(my_dens);
        plot_ste{k}(n) = std(my_dens)/sqrt(length(my_dens)-1);
        plot_normfact(n) = mean(total_vox_in_hull{k}(bin_num==n));
        prob_nonzero(n) = mean(my_dens>0);
    end
    
%     errorbar(main_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_data,plot_ste, 'LineWidth', 2);
    errorbar(main_ax, bins,plot_data{k}*100,plot_ste{k}*100, 'LineWidth', 2, 'Color', colmap(k,:));
%     plot(off_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_normfact,
%     'LineWidth', 2);
    
%     plot(main_ax, bin_size*(1:9)-bin_size/2,plot_data(1:9), 'LineWidth', 2, 'Color', c(k,:)); 
%     scatter(main_ax, bin_num*bin_size - bin_size/2 + k, density, '*', 'MarkerEdgeColor', c(k,:));
    
end

prep_figure(gcf, gca, 'legend', types, 'xlabel', 'Distance from soma (microns)', 'ylabel', 'BC contact per SAC area (um^2/um^2)');
%     
% 
set(gca, 'YTick', [0:.5:10]);


% end
% title('Bipolar to off-SAC contact density')
% legend(types)

% M = zeros(numel(D),6);
% k = 0;
% for b = 1:length(all_cell_nums)
%     for s = 1:length(dsgc_nums)
%         k = k+1;
%         M(k,1) = dsgc_nums(s);
%         M(k,2) = all_cell_nums(b);
%         M(k,3) = all_types(b);        
%         M(k,4) = D(s,b);
%         M(k,5) = H(s,b);
%         M(k,6) = X(s,b);
%     end
% end
% 
% csvwrite('fig4d.txt', M);
% 
%         




