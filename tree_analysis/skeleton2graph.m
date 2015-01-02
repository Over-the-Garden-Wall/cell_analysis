function G = skeleton2graph(skele)
    
    edge_weight = sqrt(sum((skele.nodes(skele.edges(:,1),:)-skele.nodes(skele.edges(:,2),:)).^2,2));
    
    G = sparse([skele.edges(:,1); skele.edges(:,2)], [skele.edges(:,2); skele.edges(:,1)], [edge_weight; edge_weight], size(skele.nodes,1), size(skele.nodes,1));
    
end
