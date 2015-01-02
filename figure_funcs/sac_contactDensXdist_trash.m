C = get_constants;

% types = {'t1', 't3a'};
types = {'t1', 't2', 't3a', 't3b', 't4'};

sac_nums = C.type.sure_off_sac;
% sac_nums = C.type.off_sac;
num_sacs = length(sac_nums);
min_thresh = 10000;
bin_size = 15;

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
        
    
    [total_contact{k}, total_vox_in_hull{k}] = get_contact_density_whulls(sac_nums, cell_nums);
    
    total_vox_in_hull{k} = total_vox_in_hull{k}(:)*prod(C.res)/10^9;
    cell_dist{k} = cell_dist{k}(:);
    
    
    
    bin_num = ceil(cell_dist{k}/bin_size);
    
    plot_data{k} = zeros(max(bin_num),1);
    plot_ste{k} =  zeros(max(bin_num),1);    
    prob_nonzero =  zeros(max(bin_num),1);
    plot_normfact =  zeros(max(bin_num),1);
    
    for n = 1:max(bin_num)
        tot_vol = total_vox_in_hull{k}(bin_num==n);
        plot_data{k}(n) = sum(tot_vol);
        
    end
    
%     errorbar(main_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_data,plot_ste, 'LineWidth', 2);
%     errorbar(main_ax, bin_size*(1:9)-bin_size/2,plot_data{k}(1:9),plot_ste{k}(1:9), 'LineWidth', 2, 'Color', C.colormap(k,:));
  plot(main_ax, bin_size*(1:9)-bin_size/2,plot_data{k}(1:9), 'LineWidth', 2, 'Color', C.colormap(k,:));
%     plot(off_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_normfact,
%     'LineWidth', 2);
    
%     plot(main_ax, bin_size*(1:9)-bin_size/2,plot_data(1:9), 'LineWidth', 2, 'Color', c(k,:)); 
%     scatter(main_ax, bin_num*bin_size - bin_size/2 + k, density, '*', 'MarkerEdgeColor', c(k,:));
    
end
% title('Bipolar to off-SAC contact density')
% legend(types)

% legend({'type 1', 'type 3a'})
legend({'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'})
ylabel('% contact with SAC voxels', 'FontSize', 20)
xlabel('Distance from soma (microns)', 'FontSize', 20)
set(gcf, 'Position', [0 0 640 640]);
set(gca, 'FontSize', 20);

set(gca, 'Position', [.1 .1 .8 .8]);




