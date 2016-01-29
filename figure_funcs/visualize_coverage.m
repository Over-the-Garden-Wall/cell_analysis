function visualize_coverage(cell_nums, resolution)

        num_cells = length(cell_nums);
    hull_points = cell(num_cells,2);

    min_points = [Inf Inf];
    max_points = -[Inf Inf];
    
    for k = 1:num_cells
        c_d = cell_data(cell_nums(k));
        [hull_points{k,:}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
    end
    
    for k = 1:num_cells
        for d = 1:2
            min_points(d) = min([min_points(d); hull_points{k,d}]);
            max_points(d) = max([max_points(d); hull_points{k,d}]);
        end
    end    
    
    for k = 1:num_cells
        for d = 1:2
            hull_points{k,d} = (hull_points{k,d} - min_points(d)) / resolution(d);
        end
    end
       
    im_size = ceil((max_points - min_points)./resolution);
    [im_map_y im_map_x] = meshgrid(1:im_size(2), 1:im_size(1));
    im = zeros(im_size);
    
    for k = 1:num_cells
        is_in = inpolygon(im_map_x(:), im_map_y(:), hull_points{k,1}, hull_points{k,2});
        im(is_in) = im(is_in) + 1;
    end
    
    figure; imagesc(im); colorbar
    title(num2str(mean(im(im(:)~=0))));
end
    
    