C = get_constants;

colmap = colormap('lines');
% types = {'t1', 't3a'};
% types = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'};
types = {'BC5t', 'BC5o', 'BC5i', 'XBC', 'BC6', 'BC7', 'BC89', 'RBC'};


% sac_nums = C.type.sure_off_sac;
sac_nums = C.type.on_sac;

num_sacs = length(sac_nums);
min_thresh = 10000;
bin_size = 15;
bins = bin_size/2:bin_size:bin_size*8.5;



c = colormap('Lines');


density_sum = zeros(length(bins),length(types));

% for k = 1:length(types);
%     
%     cell_nums = C.type.(types{k});
%     num_cells = length(cell_nums);
%     cell_dist = zeros(num_sacs,num_cells);    
%     for s = 1:num_sacs
%         s_dat = cell_data(sac_nums(s));
%         soma_loc = s_dat.get_midpoint(true);
%         for h = 1:num_cells
%             h_dat = cell_data(cell_nums(h));
%             mid_loc = h_dat.get_midpoint(false);
%             
%             for d = 2:3
%                 cell_dist(s,h) = cell_dist(s,h) + ...
%                     (soma_loc(d)-mid_loc(d))^2;
%             end
%             cell_dist(s,h) = sqrt(cell_dist(s,h))/1000;
%         end
%     end
%         
%     
%     [total_contact, total_vox_in_hull] = get_contact_density_whulls(sac_nums, cell_nums);
%     
%     is_valid = total_vox_in_hull > min_thresh;
%     
%     total_contact = total_contact(is_valid(:));
%     total_vox_in_hull = total_vox_in_hull(is_valid(:));
%     cell_dist = cell_dist(is_valid(:));
%     
%     
%     
%     bin_num = ceil(cell_dist/bin_size);
%     density = total_contact./total_vox_in_hull;
%     
%     
%     for n = 1:length(bins)
%         my_dens = density(bin_num==n);
%         density_sum(n,k) = mean(my_dens);
%     end    
%     
% end




% relevant_portion = [10 60];
relevant_portion = [40 80];

relevant_x = relevant_portion(1):relevant_portion(2);

strats = zeros(length(relevant_x),length(types));

[quartile_data bins full_x full_data] = get_stratification_by_dist(sac_nums, [], bins, [.25 .5 .75], true, false, 1);


mean_hull_area = zeros(length(types),1);

for k = 1:length(types);
    cell_nums = C.type.(types{k});
    for n = 1:length(cell_nums)
        cell_dat = cell_data(cell_nums(n));
        
        mean_hull_area(k) = mean_hull_area(k) + poly_area(cell_dat.hull_2d);
        
        p = cell_dat.get_surface;
        p = C.f(p(:,1));
        p = p(p>=relevant_portion(1) & p<= relevant_portion(2));
        strats(:,k) = strats(:,k) + hist(p,relevant_x)';
    end
    strats(:,k) = strats(:,k)/length(cell_nums);
    mean_hull_area(k) = mean_hull_area(k)/length(cell_nums);
    strat_norm(k) = sum(strats(:,k));
    strats(:,k) = strats(:,k)/strat_norm(k);
%     plot(strats(:,k), relevant_x, 'LineWidth', 2);
end



t_density = 1./mean_hull_area*16.5*25;
% t_density = [2233 3212 1866 3254 3005];

is_valid = full_x >= relevant_x(1) & full_x <= relevant_x(end);
full_data = full_data(is_valid,:)/num_sacs;
overlap = zeros(length(bins),length(types));

if size(full_data,1) > size(strats,1)    
    full_data = full_data(1:size(strats,1),:);
elseif size(full_data,1) < size(strats,1)
    strats = strats(1:size(full_data,1),:);
end

for k = 1:length(types)
    for n = 1:length(bins);
        
        p_s_giv_xyz = full_data(:,n); 
        p_bk_giv_z = strats(:,k)*strat_norm(k)*t_density(k); %
        
        bin_normd = sum(full_data(:,n));
        overlap(n,k) = sum(p_bk_giv_z.*p_s_giv_xyz)/sum(p_s_giv_xyz);
        
        %overlap is p_bk_giv_sxy
        
%         overlap(n,k) = .6*overlap(n,k) + .4*(2*overlap(n,k) - overlap(n,k).^2);
        % 60% of s voxels have 1 edge, 40% have 2 (according to cylindrical model)
        % basically same as 1.4x unless p is large
        
        
        overlap(n,k) = overlap(n,k) * 1.4;
%         sum(strats(:,k).*full_data(:,n))*strat_norm(k)*t_density(k);
    end
end

figure; hold all
c = colormap('Lines');
sc_h = zeros(size(overlap));
for n = 1:length(bins);
    whiteness = n/length(bins)/2;
    for k = 1:length(types);
        my_c = c(k,:)*(1-whiteness) + [1 1 1] * whiteness;
        sc_h(n,k) = scatter(overlap(n,k), density_sum(n,k), '*', 'MarkerEdgeColor', my_c);
    end
end
        
legend(sc_h(1,:), types);    



        

figure; hold on
% C = get_constants;

for k = 1:length(types)
    plot(bins, overlap(:,k)*100, 'LineWidth', 2, 'Color', C.colormap(types{k}));
end


t = get(gca, 'YLim');
% set(gca,'YLim',[0 t(2)])
set(gca,'YLim',[0 3])
    set(gca, 'TickDir', 'in');
    
    set(gca, 'FontSize', 20);
    
    set(gca, 'Box', 'off');
    
%     set(gca, 'YTick', [.5 1 1.5 2]);
    
    
    set(gcf, 'Position', [0 0 640 640]);
    set(gca, 'Position', [.1 .1 .8 .8]);
    
    


