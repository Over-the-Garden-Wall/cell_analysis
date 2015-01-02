function [overlap, weighted_overlap, total_area] = find_stratification_overlap(bip_cells, target_cell, Q)
    
    res = [16.5 16.5 25];
    
    f = inline(' (30-108)/(2.7-1.5)/10^4.*(x-2.7*10^4)+30');
                
    
    %get stratification for target (jamb) cell
    
    fn = ['./surface_points/cell_' num2str(target_cell) '_surface.mat'];        
    
%         if exist(fn,'file')
    load(fn);
    for k = 1:3
        surface_points(:,k) = surface_points(:,k)*res(k);
    end                

    surface_points = surface_points*Q';
    
%     [min_z min_ind] = min(surface_points(:,1));
%     prox_point = mean(surface_points(surface_points(:,1)==min_z,:),1);
%     dist_from_prox = (surface_points(:,1) - prox_point(1)).^2 + ...
%         (surface_points(:,2) - prox_point(2)).^2 + ...
%         (surface_points(:,3) - prox_point(3)).^2;
% 
%     [max_dist max_ind] = max(dist_from_prox);
%     dist_point = surface_points(max_ind,:);
%     dist_vec = dist_point - prox_point;
%     dist_vec = dist_vec(2:3) / sqrt(sum(dist_vec(2:3).^2));
%     
%     dist = (surface_points(:,2) - prox_point(2)).*dist_vec(1) + ...
%         (surface_points(:,3) - prox_point(3)).*dist_vec(2);
%     
%     depth = f(surface_points(:,1));
%     depth_ticks = min(depth):max(depth);
%     depth = floor(depth-min(depth)+1);
%     
%     my_plot = zeros(max(depth),100);
%         
%     for n = 1:length(depth)
%         my_plot(depth(n),dist(n)) = my_plot(depth(n),dist(n))+1;
%     end
%     
%     figure;imagesc(my_plot); colorbar;
%     title(['number of surface points by depth and distance for ' num2str(10010)])

    min_point = min(surface_points);
    max_point = max(surface_points);
    
    num_bins = [100 100 100];
    binned_points = zeros(num_bins);
    
    for n = 1:size(surface_points,1);
        point_bin = floor((surface_points(n,:) - min_point)./(max_point-min_point).*(num_bins-1)+1);
        binned_points(point_bin(1), point_bin(2), point_bin(3)) = 1 + ...
            binned_points(point_bin(1), point_bin(2), point_bin(3));
    end
        
        
    
    num_cells = length(bip_cells);
    
    overlap = zeros(num_cells,1);
    weighted_overlap = zeros(num_cells,1);
    total_area = zeros(num_cells,1);
        
    for n = 1:num_cells
        
        fn = ['./surface_points/cell_' num2str(bip_cells(n)) '_surface.mat'];
        if exist(fn,'file')
            load(fn);
            for k = 1:3
                surface_points(:,k) = surface_points(:,k)*res(k);
            end                
            num_points = size(surface_points,1);
            disp(['loaded ' num2str(num_points) ' points from cell ' num2str(bip_cells(n))])
            
            surface_points = surface_points*Q';
            
            total_area(n) = num_points;
            
            point_bin = zeros(size(surface_points));
            is_valid = true(num_points,1);
            for k = 1:3
                point_bin(:,k) = floor((surface_points(:,k) - min_point(k))/(max_point(k)-min_point(k))*(num_bins(k)-1)+1);
                is_valid = is_valid & point_bin(:,k) > 0 & point_bin(:,k) <= num_bins(k);
            end
            
            bin_val = zeros(num_points,1);
            for k = find(is_valid')                
                bin_val(k) = binned_points(point_bin(k,1), point_bin(k,2), point_bin(k,3));
            end
                        
            overlap(n) = sum(bin_val > 0);
            weighted_overlap(n) = sum(bin_val);
            
            
        end
        
    end
    
end
    
    
    
    
    