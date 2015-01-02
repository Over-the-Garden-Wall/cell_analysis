function [density_matrix, overlap_matrix, cn1_density, cn2_density] = connectivity_histogram(cn1, cn2, min_overlap)

    nc1 = length(cn1);
    nc2 = length(cn2);

    overlap_matrix = zeros(nc1, nc2);
    connectivity_matrix = zeros(nc1, nc2);
    
    
    for m = 1:nc1
        c_d1 = cell_data(cn1(m));
        my_conns = double(c_d1.contacts);
        for n = 1:nc2
            c_d2 = cell_data(cn2(n));
            
            h1 = c_d1.hull_2d;
            [h1(:,1) h1(:,2)] = poly2cw(h1(:,1), h1(:,2));
            
            h2 = c_d2.hull_2d;
            [h2(:,1) h2(:,2)] = poly2cw(h2(:,1), h2(:,2));
            
            h = [];
            try
                [h(:,1), h(:,2)] = polybool('intersection', h1(:,1), h1(:,2), h2(:,1), h2(:,2));
                ol = polyarea(h(:,1), h(:,2));
            catch ME
                ol = 0;
            end
            
            overlap_matrix(m,n) = ol;
            
            conn_is_c2 = my_conns(1,:) == cn2(n);
            
            connectivity_matrix(m,n) = sum(my_conns(2,conn_is_c2));
        end
    end
    
    overlap_matrix(overlap_matrix < min_overlap) = 0;
    
    density_matrix = connectivity_matrix ./ overlap_matrix;
    density_matrix(isinf(density_matrix)) = NaN;
    
    cn2_density = sum(connectivity_matrix)./sum(overlap_matrix);
    cn1_density = sum(connectivity_matrix,2)./sum(overlap_matrix,2);
    
    figure; hist(cn1_density(~isnan(cn1_density) & ~isinf(cn1_density)),20);
    figure; hist(cn2_density(~isnan(cn2_density) & ~isinf(cn2_density)),20);
    
end
    