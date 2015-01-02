function p = tree2points(nodes, edges, points_per_unit_length)
    
    total_length = get_skele_length(nodes,edges);
    
    p = zeros(floor(total_length/points_per_unit_length),3);
    
    remainder = 0;
    k = 0;
    for n = 1:length(edges)
        edge_length = sqrt(sum((nodes(edges(n,1),:)-nodes(edges(n,2),:)).^2));
        num_ps = floor((edge_length+remainder)/points_per_unit_length);
        for m = 1:num_ps
            alpha = (m-remainder)*points_per_unit_length/edge_length;
            p(k+m,:) = nodes(edges(n,1),:)*(1-alpha) + nodes(edges(n,2),:)*alpha;
        end
        k = k+num_ps;
        
        if k >  floor(total_length/points_per_unit_length)
            k=k;
        end
        
        remainder = edge_length - num_ps*points_per_unit_length + remainder;
    end
    
end