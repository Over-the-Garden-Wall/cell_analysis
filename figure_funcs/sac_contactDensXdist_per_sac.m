C = get_constants;

% types = {'t1', 't3a'};
% types = {'t1', 't2', 'sure_t3a', 't3b', 'sure_t4'};
types = {'t1', 't2', 't3a', 't3b', 't4'};
num_types = length(types);
% types = {'t3a'};
sac_nums = C.type.sure_off_sac;
% sac_nums = C.type.off_sac;
num_sacs = length(sac_nums);
min_thresh = 10000;
bins = C.sac_bins;
bin_size = bins(2)-bins(1);
num_bins = length(bins);

% figure; off_ax = gca; hold all

c = colormap('Lines');


hulls = cell(num_types,1);
for k = 1:length(types);
    cell_nums = C.type.(types{k});
    num_cells = length(cell_nums);
    hulls{k} = cell(num_cells,1);
    
    for c = 1:num_cells
        cell_dat = cell_data(cell_nums(c));
        hulls{k}{c} = cell_dat.hull_2d;
    end
end

cell_dist = cell(num_sacs,num_types);
bin_num = cell(num_sacs,num_types);
cell_contact = cell(num_sacs,num_types);
bin_hulls = cell(num_sacs, num_types, num_bins);
points_in_bin = zeros(num_sacs, num_types, num_bins);    

total_dens = zeros(num_sacs, num_types, num_bins);

for s = 1:num_sacs
    s_dat = cell_data(sac_nums(s));
    soma_loc = s_dat.get_midpoint(true);
    p = s_dat.get_surface;
    p = p(:,2:3);
    
    z = double(s_dat.contacts(3,:));
    is_valid = C.f(z) > 10 & C.f(z) < 70;
    
    for k = 1:length(types);
    
        cell_nums = C.type.(types{k});
        num_cells = length(cell_nums);
        cell_dist{s,k} = zeros(num_cells,1);    
        cell_contact{s,k} = zeros(num_cells,1);    
        
        for h = 1:num_cells
            h_dat = cell_data(cell_nums(h));
            mid_loc = h_dat.get_midpoint(false);
            
            cell_dist{s,k}(h) = sqrt((soma_loc(2)-mid_loc(2))^2 + (soma_loc(3)-mid_loc(3))^2)/1000;
            
            for d = 2:3
                cell_dist{s,k}(h) = cell_dist{s,k}(h) + ...
                    (soma_loc(d)-mid_loc(d))^2;
            end
            cell_dist{s,k}(h) = sqrt(cell_dist{s,k}(h))/1000;
  
            is_me = double(s_dat.contacts(1,:)) == cell_nums(h);
            
            cell_contact{s,k}(h) = sum(double(s_dat.contacts(2,is_me & is_valid)));
%             
%             if s_dat.contact_map.isKey(cell_nums(h))
%                 cell_contact{s,k}(h) = s_dat.contact_area(s_dat.contact_map(cell_nums(h)));
%             else
%                 cell_contact{s,k}(h) = 0;
%             end                                        
        
        end
        
        bin_num{s,k} = ceil((cell_dist{s,k}-bins(1))/bin_size+.5);
  
        for n = 1:num_bins
            cells_in_bin = find(bin_num{s,k}==n);
            for c = cells_in_bin'
                if isempty(bin_hulls{s,k,n})
                    bin_hulls{s,k,n} = hulls{k}{c};
                else
                    temphull = [];
                    [temphull(:,1), temphull(:,2)] = polybool('union', bin_hulls{s,k,n}(:,1), bin_hulls{s,k,n}(:,2), hulls{k}{c}(:,1), hulls{k}{c}(:,2));
                    bin_hulls{s,k,n} = temphull;
                end
            end
        end
    
        for n = 1:num_bins
            if ~any(bin_num{s,k}(:)==n)
                points_in_bin(s,k,n) = 0;            
                total_dens(s,k,n) = nan;
            else
%                 tic
                points_in_bin(s,k,n) = sum(rough_2d_inpolygon(p,bin_hulls{s,k,n},[150 150]));
%                 points_in_bin(s,k,n) = sum(inpolygon(p(:,1),p(:,2), bin_hulls{s,k,n}(:,1), bin_hulls{s,k,n}(:,2)));            
%                 toc
                total_dens(s,k,n) = sum(cell_contact{s,k}(bin_num{s,k}(:)==n))/points_in_bin(s,k,n);        
            end
        end
        
        
        
    end
    

end

total_dens(total_dens>1) = 1;

mean_density = zeros(num_types, num_bins);
ste_density = zeros(num_types, num_bins);


figure; main_ax = gca; hold all;

for k = 1:num_types;
    for n = 1:num_bins        
        is_valid = ~isnan(total_dens(:,k,n));
        mean_density(k,n) = mean(total_dens(is_valid,k,n));
        ste_density(k,n) = std(total_dens(is_valid,k,n))/sqrt(sum(is_valid)-1);
    end
    errorbar(main_ax, bins,mean_density(k,:)*100,ste_density(k,:)*100, 'LineWidth', 2, 'Color', C.colormap(k,:));    
end
        
    
%     errorbar(main_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_data,plot_ste, 'LineWidth', 2);
%     plot(off_ax, bin_size*(1:max(bin_num))-bin_size/2,plot_normfact,
%     'LineWidth', 2);
    
%     plot(main_ax, bin_size*(1:9)-bin_size/2,plot_data(1:9), 'LineWidth', 2, 'Color', c(k,:)); 
%     scatter(main_ax, bin_num*bin_size - bin_size/2 + k, density, '*', 'MarkerEdgeColor', c(k,:));
    
    
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




