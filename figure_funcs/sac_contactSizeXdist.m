C = get_constants;

% types = {'t1', 't3a'};
% types = {'t1', 't2', 'sure_t3a', 't3b', 'sure_t4'};
% types = {'t1', 't2', 't3a', 't3b', 't4'};
% types = {'t2', 't3a'};
% sac_nums = C.type.sure_off_sac;
sac_nums = C.type.on_sac;

colmap = C.colormap;
colmap = colormap('Lines');
% types = {'xbc', 't5w', 't5h', 't5l', 't6', 't7', 't89', 'tRBC'};
types = {'t5w', 't5l', 't5h', 'xbc', 't6', 't7', 't89', 'tRBC'};

% g = []; for k = 1:11; g{end+1} = ['g' num2str(k)]; end; types = g;

% sac_nums = C.type.off_sac;
num_sacs = length(sac_nums);
min_thresh = 10000;
bins = C.sac_bins;
bin_size = bins(2)-bins(1);
num_bins = length(bins);

figure; main_ax = gca; hold all;
% figure; off_ax = gca; hold all

c = colormap('Lines');


contact_total = zeros(num_bins, length(types));
contact_count = zeros(num_bins, length(types));

for k = 1:length(types);
    
    cell_nums = C.type.(types{k});
    num_cells = length(cell_nums);
%     cell_dist{k} = zeros(num_sacs,num_cells);    
    for s = 1:num_sacs
        s_dat = cell_data(sac_nums(s));
        soma_loc = s_dat.get_midpoint(true);
        for h = 1:num_cells
            h_dat = cell_data(cell_nums(h));
            mid_loc = h_dat.get_midpoint(false);
            cell_dist = sqrt(sum((soma_loc(2:3)-mid_loc(2:3)).^2))/1000;
            
            [dummy, bin_num] = min(abs(bins-cell_dist));
            
            is_me = s_dat.contacts(1,:) == cell_nums(h);
            contact_count(bin_num, k) = contact_count(bin_num, k) + sum(is_me);
            contact_total(bin_num, k) = contact_total(bin_num, k) + sum(s_dat.contacts(2,is_me));
            
        end
    end
        
    
    plot(main_ax, bins, contact_total(:,k)./contact_count(:,k), 'LineWidth', 2, 'Color', c(k,:)); 
%     scatter(main_ax, bin_num*bin_size - bin_size/2 + k, density, '*', 'MarkerEdgeColor', c(k,:));
    
end
% title('Bipolar to off-SAC contact density')
% % legend(types)
% 
% M = zeros(numel(D),6);
% k = 0;
% for b = 1:length(all_cell_nums)
%     for s = 1:length(sac_nums)
%         k = k+1;
%         M(k,1) = sac_nums(s);
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
%     
% % legend({'type 1', 'type 3a'})
% legend({'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'})
% ylabel('BC contact per SAC area (um^2/um^2)', 'FontSize', 20)
% xlabel('Distance from soma (microns)', 'FontSize', 20)
% 
% set(gca, 'YTick', [0:.5:10]);

set(gcf, 'Position', [0 0 640 640]);
set(gca, 'FontSize', 20);

set(gca, 'Position', [.1 .1 .8 .8]);




