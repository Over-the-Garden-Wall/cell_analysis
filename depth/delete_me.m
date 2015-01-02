function [depthXdist soma_points] = count_depth_by_distance(cell_nums, max_dist, given_soma)
    
    num_cells = length(cell_nums);
    soma_points = zeros(num_cells,3);
    
    f = @z2ipldepth;
    
    depthXdist = zeros(100, max_dist);
    
    for n = 1:num_cells
        fn = ['./surface_points_trans/cell_' num2str(cell_nums(n)) '_surface.mat'];
        if exist(fn,'file')
            load(fn);
            
            if ~exist('given_soma','var') || isempty(given_soma)
                [min_point, min_ind] = max(surface_points(:,1));
                warning('soma detection only works for cells with cell body in the INL')
                soma_points(n,:) = surface_points(min_ind,:);
            else
                soma_points(n,:) = surface_points(min_ind,:);
            end
            
            num_points = size(surface_points,1);
            dist = sqrt(sum((surface_points - ones(num_points,1)*soma_points(n,:)).^2,2));

            
            dist = ceil(dist/1000); %distance in microns
            depth = ceil(f(surface_points(:,1)));
            
            is_valid = depth > 0 & depth < 100 & dist <= max_dist & dist > 0;
            dist = dist(is_valid);
            depth = depth(is_valid);
            
            for k = 1:length(depth)
                depthXdist(depth(k),dist(k))=depthXdist(depth(k),dist(k))+1;
                
            end
        end
    end        
end