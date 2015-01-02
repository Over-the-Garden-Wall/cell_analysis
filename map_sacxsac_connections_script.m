function [sac_conn_map diff_map dist_map soma_dist] = map_sacxsac_connections_script
    C = get_constants;
    max_dist = 200;
    
    load(C.conn_loc);
    
    conns = pick_multi_conns(conns, {C.type.off_sac, C.type.off_sac}, true);
    
    num_cells = length(C.type.off_sac);
    somata = zeros(num_cells,3);
    for k = 1:num_cells
        somata(k,:) = get_soma_loc(C.type.off_sac(k));
    end
    
    sac_conn_map = zeros(180, max_dist, max_dist);
    diff_map = zeros(180, max_dist);
    dist_map = zeros(180, max_dist*2);
    soma_dist = zeros(num_cells);
    
    
    
    for k = 1:num_cells
        for l = 1:num_cells
            hyp = sqrt((somata(k,2) - somata(l,2))^2 + (somata(k,3) - somata(l,3))^2);
            soma_dist(k,l) = hyp;
            
            conns{k,l} = double(conns{k,l});
            for n = 1:size(conns{k,l},2);
                dist1 = sqrt((conns{k,l}(5,n) - somata(k,2))^2 + (conns{k,l}(6,n) - somata(k,3))^2);
                dist2 = sqrt((conns{k,l}(5,n) - somata(l,2))^2 + (conns{k,l}(6,n) - somata(l,3))^2);
                
                gamma = acos((dist1^2 + dist2^2 - hyp^2)/(2*dist1*dist2));
                gamma = ceil(gamma/pi*180);
                if gamma == 0;
                    gamma = 1;
                end
                dist1 = ceil(dist1/1000);
                dist2 = ceil(dist2/1000);                
                
                dist_map(gamma, ceil(hyp/1000)) = dist_map(gamma, ceil(hyp/1000)) + 1;
                sac_conn_map(gamma, dist1, dist2) = sac_conn_map(gamma, dist1, dist2)+conns{k,l}(3,n);
            end
        end
    end
    
    for k = 1:max_dist
        for l = k+1:max_dist
            d = l-k+1;
            
            diff_map(:,d) = diff_map(:,d) + sac_conn_map(:,k,l);
        end
    end
    
            
            
    
    
    
end