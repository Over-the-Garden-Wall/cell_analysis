
C = get_constants;


for sac_nums = C.type.on_sac


% types = {'t1', 't3a'};
% types = {'t1', 't2', 'sure_t3a', 't3b', 'sure_t4'};
% types = {'t1', 't2', 't3a', 't3b', 't4'};
% types = {'t2', 't3a'};
% sac_nums = C.type.sure_off_sac;
% sac_nums = C.type.on_sac;

colmap = C.colormap;
colmap = colormap('Lines');
% types = {'xbc', 't5w', 't5h', 't5l', 't6', 't7', 't89', 'tRBC'};
types = {'BC5t', 'BC5i', 'BC7'};
% types = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'};
% types = {'BC2', 'BC3a', 'AC14'};

% g = []; for k = 1:11; g{end+1} = ['g' num2str(k)]; end; types = g;

% sac_nums = C.type.off_sac;
num_sacs = length(sac_nums);
min_thresh = 10000;
bins = C.sac_bins;
bin_size = bins(2)-bins(1);
num_bins = length(bins);

total_mean = zeros(length(types),1);
total_ste = zeros(length(types),1);

% figure; off_ax = gca; hold all

c = colormap('Lines');

mof = ones(1,100);
for k = 1:length(types)
    mof(k) = estimate_mosaic_overlap(C.type.(types{k}));
end
mof(9) = 1; 

plot_data = [];
plot_ste = [];
total_contact = [];
total_vox_in_hull = [];
cell_dist = [];

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
        
    
    [total_contact{k}, total_vox_in_hull{k}] = get_contact_density_whulls(sac_nums, cell_nums, [], []);
    
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
    total_mean(k) = mean(density(:));
    total_ste(k) = std(density(:))/sqrt(length(density(:))-1);
    
    if ~isnan(plot_data{1}(1))
    
        if k ==1
            figure; main_ax = gca; hold all;
        end
%     errorbar(main_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_data,plot_ste, 'LineWidth', 2);
    errorbar(main_ax, bins,plot_data{k}*100,plot_ste{k}*0, 'LineWidth', 2, 'Color', C.colormap(types{k}));
%     plot(off_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_normfact,
%     'LineWidth', 2);
    end
%     plot(main_ax, bin_size*(1:9)-bin_size/2,plot_data(1:9), 'LineWidth', 2, 'Color', c(k,:)); 
%     scatter(main_ax, bin_num*bin_size - bin_size/2 + k, density, '*', 'MarkerEdgeColor', c(k,:));
    
end
% title('Bipolar to off-SAC contact density')
% legend(types)

M = zeros(numel(D),6);
k = 0;
for b = 1:length(all_cell_nums)
    for s = 1:length(sac_nums)
        k = k+1;
        M(k,1) = sac_nums(s);
        M(k,2) = all_cell_nums(b);
        M(k,3) = all_types(b);        
        M(k,4) = D(s,b);
        M(k,5) = H(s,b);
        M(k,6) = X(s,b);
    end
end

% csvwrite('fig4d.txt', M);

        
prep_figure(gcf, gca, 'legend', types, 'xlabel', 'Distance from soma (microns)', 'ylabel', 'BC contact per SAC area (um^2/um^2)');
    

set(gca, 'YLim', [0 6], 'YTick', [0:.5:10]);
set(gca, 'XLim', [bins(1) bins(end)], 'XTick', bins);

title(num2str(sac_nums));

end

% 
% figure;
% barplot_werror(total_mean, total_ste);
% prep_figure(gcf,gca);

