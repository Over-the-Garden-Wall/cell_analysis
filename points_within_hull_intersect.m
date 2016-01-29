function olp = points_within_hull_intersect(src_ids, dest_ids)
    
    P_SPARSITY = 100;

    all_ids = {src_ids, dest_ids};
    num_types = [length(src_ids), length(dest_ids)];
            
    type_hulls = cell(1, 2);
    
    h = cell(1,2);
    for t = 1:2
        type_hulls{t} = cell(num_types(t),2);
        for k = 1:num_types(t);
            for n = all_ids{t}{k}
                c_d = cell_data(n);
                [h{1}, h{2}] = poly2cw(c_d.hull_2d(:,1),c_d.hull_2d(:,2));
                [type_hulls{t}{k,:}] = polybool('union', h{1}, h{2}, type_hulls{t}{k,1}, type_hulls{t}{k,2});
            end
        end
    end
    
    olp = zeros(num_types);
    roi = cell(1,2);
    for k = 1:num_types(1)
        for kk = 1:num_types(2)
            [roi{:}] = polybool('intersection', ...
                type_hulls{1}{k,1}, type_hulls{1}{k,2}, ...
                type_hulls{2}{kk,1}, type_hulls{2}{kk,2});

            for c = src_ids{k}
                c_d = cell_data(c);
                p = c_d.get_surface;
                p = p(round(P_SPARSITY/2):P_SPARSITY:end,2:3);
                olp(k,kk) = olp(k,kk) + sum(inpolygon(p(:,1), p(:,2), roi{1}, roi{2}));
            end            
        end        
    end
    
    olp = olp*P_SPARSITY;
    
end