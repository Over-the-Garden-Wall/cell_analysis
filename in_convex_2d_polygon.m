function p = in_convex_2d_polygon(p, poly)

    num_edges = size(poly,1)-1;
    poly_center = mean(poly);
    
    for n = 1:num_edges
        p_c = [p(:,1) - poly(n,1), p(:,2) - poly(n,2)];
        poly_c = [poly(:,1) - poly(n,1), poly(:,2) - poly(n,2)];
        center_c = [poly_center(1) - poly(n,1), poly_center(2) - poly(n,2)];
        
        v = [poly_c(n+1,2), -poly_c(n+1,1)];
        
        p_d = p_c(:,1)*v(1) + p_c(:,2)*v(2);
        
        good_sign = center_c(1)*v(1) + center_c(2)*v(2);
        p = p(p_d*good_sign >= 0, :);
        p_c = p_c(p_d*good_sign >= 0, :);
    end
        
end