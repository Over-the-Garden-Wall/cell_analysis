function [cell_plots plot_stds] = depth_by_distance_from_soma(cell_nums, max_dist)
    
    soma_threshold = 200;

    num_cells = length(cell_nums);

    cell_plots = zeros(num_cells,max_dist);
    plot_stds = zeros(num_cells,max_dist);
    is_valid = true(num_cells,1);
    f = @z2ipldepth;        
    
    for n = 1:num_cells
        fn = ['./surface_points_trans/cell_' num2str(cell_nums(n)) '_surface.mat'];
        if exist(fn,'file')
            load(fn);
            [min_point] = max(surface_points(:,1));
            
            soma_point = mean(surface_points(surface_points(:,1)>= min_point - soma_threshold, :));

            num_points = size(surface_points,1);
            dist = sqrt(sum((surface_points(:,2:3) - ones(num_points,1)*soma_point(2:3)).^2,2));

            
            dist = ceil(dist/1000);
            
            for k = 1:min(max(dist),max_dist)
                cell_plots(n,k) = mean(f(surface_points(dist == k,1)));
                plot_stds(n,k) = std(f(surface_points(dist == k,1)));
            end
        else
            is_valid(n) = false;            
        end
    end        
    
    cell_plots = cell_plots(is_valid,:);
    plot_stds(is_valid,:);
end