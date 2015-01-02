function full_len = get_skele_length(nodes,edges)

    
    full_len = 0;
    for n = 1:size(edges,1)
        
        full_len = full_len + sqrt(sum((nodes(edges(n,1),:) - nodes(edges(n,2),:)).^2));
%         
%         plot(xyaxis, coords(id2index([src_edge dest_edge]),1), coords(id2index([src_edge dest_edge]),2));
%         plot(yzaxis, coords(id2index([src_edge dest_edge]),2), coords(id2index([src_edge dest_edge]),3));
        
        
    end
end
