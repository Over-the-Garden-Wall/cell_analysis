function total_length = calc_skele_length(nodes, edges)
    res = [16.5 16.5 25]*2;
    
    for d = 1:3;
        nodes(:,d) = nodes(:,d)*res(d);
    end
    
    total_length = 0;
    for n = 1:size(edges,1);
        p1 = nodes(edges(n,1),:);
        p2 = nodes(edges(n,2),:);
        
        edge_dist = sqrt(sum((p1-p2).^2));
        
        total_length = total_length+edge_dist;
    end
end
        