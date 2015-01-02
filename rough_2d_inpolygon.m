function is_in = rough_2d_inpolygon(p, poly, grid_size)
    %fast estimation of inpolygon using a coarse grid O(n) + O(prod(grid_size))
    
%     warning('not debugged');
    
    poly_min = min(poly);
    poly_max = max(poly)-poly_min;
    
    
    for k = 1:size(p,2)
        p(:,k) = (p(:,k)-poly_min(k))/poly_max(k);
        poly(:,k) = (poly(:,k)-poly_min(k))/poly_max(k);
    end
    
    step_size = 1./grid_size;
    [X Y] = meshgrid(step_size(1):step_size(1):1, step_size(2):step_size(2):1);
    X = X-step_size(1)/2;
    Y = Y-step_size(2)/2;
    
    cube_is_in = inpolygon(X(:),Y(:),poly(:,1), poly(:,2));
    cube_is_in(end+1) = false;
        
    x_bin = ceil(p(:,1)*grid_size(1));
    y_bin = ceil(p(:,2)*grid_size(2));
    
    out_of_bounds = x_bin <= 0 | x_bin > grid_size(1) | y_bin <= 0 | y_bin > grid_size(2);
    
    point_membership = (ceil(p(:,1)*grid_size(1))-1)*grid_size(2) + ceil(p(:,2)*grid_size(2));
    point_membership(out_of_bounds) = length(cube_is_in);
    
    is_in = cube_is_in(point_membership);
end
    