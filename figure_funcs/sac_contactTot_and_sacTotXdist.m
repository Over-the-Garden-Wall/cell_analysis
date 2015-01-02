C = get_constants;

% types = {'t1', 't3a'};
% types = {'t1', 't2', 'sure_t3a', 't3b', 'sure_t4'};
types = {'t1', 't2', 't3a', 't3b', 't4'};
% types = {'t2', 't3a'};
sac_nums = C.type.sure_off_sac;
% sac_nums = C.type.off_sac;
num_sacs = length(sac_nums);
min_thresh = 10000;
bins = C.sac_bins;
bin_size = bins(2)-bins(1);
num_bins = length(bins);

figure; main_ax = gca; hold all;
% figure; off_ax = gca; hold all

c = colormap('Lines');

for k = 1:length(types);
    
    cell_nums = C.type.(types{k});
    num_cells = length(cell_nums);
    cell_dist{k} = zeros(num_sacs,num_cells);    
    for s = 1:num_sacs
        s_dat = cell_data(sac_nums(s));
        soma_loc = s_dat.get_midpoint(true);
        for h = 1:num_cells
            h_dat = cell_data(cell_nums(h));
            mid_loc = h_dat.get_midpoint(false);
            
            for d = 2:3
                cell_dist{k}(s,h) = cell_dist{k}(s,h) + ...
                    (soma_loc(d)-mid_loc(d))^2;
            end
            cell_dist{k}(s,h) = sqrt(cell_dist{k}(s,h))/1000;
        end
        
    end
    
    total_hull = [];
    
    
    total_contact{k} = get_contact_density_whulls(sac_nums, cell_nums);
    
%     is_valid = total_vox_in_hull{k} > min_thresh;
    total_contact{k} = total_contact{k}(:);
%     total_contact{k} = total_contact{k}(is_valid(:))*.2915;
    cell_dist{k} = cell_dist{k}(:);
    
    
%     
%     
    % sac external edges per surface voxel = 2.2488, though we're not using
    % that....
    
%      mof = [1.2839    1.2945    1.1136    1.0828    1.1160];
    mof = [1.3952    1.3325    1.1761    1.1139    1.1452];
%     mof = [1.3952    1.3325    1.9112    1.2817    1.2875];
%      mof = [1.4639    1.6033    1.7456    1.3385    1.2451];
%      mof = [1.4639    1.6033    1.5456    1.2628    1.4905];
%  mof = [1.1351    1.1268    1.0796    1.0584    1.0481];
%     mof = ones(1,5);
    
    bin_num = ceil((cell_dist{k}-bins(1))/bin_size+.5);
%    density = total_contact{k}./total_vox_in_hull{k}/.410*mof(k);
    
    cont_total{k} = zeros(num_bins,1);
    sac_total{k} =  zeros(num_bins,1);
    
    
    for n = 1:num_bins
        in_my_bin = bin_num==n;
        
        cont_total{k}(n) = sum(total_contact{k}(is_in_my_bin));
        sac_total{k} =  zeros(num_bins,1);
    
        
        my_dens = density(bin_num==n);
%         county{k}(n) = sum(bin_num==n);
        plot_data{k}(n) = mean(my_dens);
        plot_ste{k}(n) = std(my_dens)/sqrt(length(my_dens)-1);
        plot_normfact(n) = mean(total_vox_in_hull{k}(bin_num==n));
        prob_nonzero(n) = mean(my_dens>0);
    end
    
%     errorbar(main_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_data,plot_ste, 'LineWidth', 2);
    errorbar(main_ax, bins,plot_data{k}*100,plot_ste{k}*100, 'LineWidth', 2, 'Color', C.colormap(k,:));
%     plot(off_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_normfact,
%     'LineWidth', 2);
    
%     plot(main_ax, bin_size*(1:9)-bin_size/2,plot_data(1:9), 'LineWidth', 2, 'Color', c(k,:)); 
%     scatter(main_ax, bin_num*bin_size - bin_size/2 + k, density, '*', 'MarkerEdgeColor', c(k,:));
    
end
% title('Bipolar to off-SAC contact density')
% legend(types)

% legend({'type 1', 'type 3a'})
legend({'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'})
ylabel('BC contact per SAC area (um^2/um^2)', 'FontSize', 20)
xlabel('Distance from soma (microns)', 'FontSize', 20)

set(gca, 'YTick', [0:.5:10]);

set(gcf, 'Position', [0 0 640 640]);
set(gca, 'FontSize', 20);

set(gca, 'Position', [.1 .1 .8 .8]);




