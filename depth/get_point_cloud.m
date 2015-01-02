function p = get_point_cloud(cell_nums, density)


    surface_dir = './surface_points_trans/';
    num_cells = length(cell_nums);
    
    fns = cell(num_cells,1);
    total_points = 0;
    for n = 1:num_cells
        fns{n} = [surface_dir 'cell_' num2str(cell_nums(n)) '_surface.mat'];
        fn_data = whos('-file', fns{n});
        total_points = total_points + ceil(fn_data.size(1)*density);
    end
    
    p = zeros(total_points,3);
    k = 0;
    for n = 1:num_cells
        load(fns{n});
        p_to_add = ceil(size(surface_points,1)*density);
        p(k + (1:p_to_add),:) = surface_points(floor(1:1/density:end),:);
        k = k+p_to_add;
    end
    
end