function [path_info, coverage_perc coverage_hull] = find_most_covered_dendrites(cell_nums, types, dendrite_region, threshold)

    if ~exist('dendrite_region', 'var') || isempty(dendrite_region)
        dendrite_region = [0 Inf];
    end

    C = get_constants;
    num_cells = length(cell_nums);
    
    path_info = [];
    coverage_perc = [];
    
    coverage_hull = find_covered_area(types);
    
    for k = 1:num_cells
        c = cell_nums(k);
        fn = ['./sac_dendrite_paths/path_' num2str(c) '.mat'];
        
        if exist(fn,'file')
            paths = [];
            load(fn);
            
            num_paths = length(paths);
            path_coverage = zeros(num_paths,1);
            
            for n = 1:length(paths)
                
                dist = sqrt(sum((paths{n}(1:end-1,:)-paths{n}(2:end,:)).^2,2));
                total_dist = [0; cumsum(dist)];
                
                is_valid = total_dist >= dendrite_region(1) & total_dist <= dendrite_region(2);
                
                path_coverage(n) = mean(inpolygon(paths{n}(is_valid,2), paths{n}(is_valid,3), coverage_hull(:,1), coverage_hull(:,2)));
%                 path_coverage(n) = path_coverage(n)/length(paths{n});


%                 figure; plot(coverage_hull(:,1), coverage_hull(:,2)); hold on; scatter(paths{n}(:,2), paths{n}(:,3)); title([num2str(n) ': ' num2str(path_coverage(n))]);

            end
               
            is_covered = path_coverage>=threshold;
            coverage_perc = [coverage_perc; path_coverage(is_covered)];
            path_info = [path_info; c*ones(sum(is_covered),1), find(is_covered)];
            
        end
    end

end

    