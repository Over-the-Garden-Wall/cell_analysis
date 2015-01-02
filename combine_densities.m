function density = combine_densities(cell_nums, varargin)

    

    max_size = [0 0 0];
    num_cells = length(cell_nums);
    indi_dens = cell(num_cells,1);
    
    for k = 1:num_cells
        indi_dens{k} = get_density_all(cell_nums(k), varargin{:});
            
        max_size = max([max_size; size(indi_dens{k})]);
    end
        
    density = zeros(max_size);
    for k = 1:num_cells
        density(1:size(indi_dens{k},1), 1:size(indi_dens{k},2),  1:size(indi_dens{k},3)) = ...
            density(1:size(indi_dens{k},1), 1:size(indi_dens{k},2),  1:size(indi_dens{k},3)) + ...
            indi_dens{k};
    end
end