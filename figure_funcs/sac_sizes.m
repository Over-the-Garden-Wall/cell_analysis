sac_types = {'on_sac', 'off_sac'};
min_dist_from_soma = 80000;
points_per_sac = 8;

min_dist_from_furthest = 30000;

C = get_constants;


cell_dist = cell(length(sac_types),1);
cell_max_dist = cell(length(sac_types),1);

for t = 1:2
    
    cns = C.type.(sac_types{t});
    
    
    for cn = 1:length(cns)
        c = cns(cn);
        c_d = cell_data(c);
        soma_loc = c_d.get_midpoint(true);
        
        p = c_d.get_surface;
        d = sqrt((p(:,2)-soma_loc(2)).^2 + (p(:,3)-soma_loc(3)).^2);
        p(d < min_dist_from_soma,:) = [];
        
        
        d_p = zeros(points_per_sac,1);
        for k = 1:points_per_sac;
            d = sqrt((p(:,2)-soma_loc(2)).^2 + (p(:,3)-soma_loc(3)).^2);
        
            [max_val, max_ind] = max(d);
        
            d_p(k) = max_val;
            
            d = sqrt((p(:,2)-p(max_ind,2)).^2 + (p(:,3)-p(max_ind,3)).^2);
            p(d < min_dist_from_furthest,:) = [];
            if isempty(p)
                break
            end
        end
        
        cell_dist{t}(cn) = mean(d_p(1:k));
        cell_max_dist{t}(cn) = max(d_p(1:k));
    end
            
end
        
        