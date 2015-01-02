function EL = graph2edgeList(G)

    parents = find(any(G),1,'first');
    
    EL = [];
    
    while any(G(:)); 
        if isempty(parents)
            old_parents = find(any(G),1,'first');
        else        
            old_parents = parents;
        end
        parents = [];
        
        for p = old_parents
            new_parents = find(G(p,:));
            for k = 1:length(new_parents)
                EL = [EL; p, new_parents(k)];
            end
            parents = [parents new_parents];            
        end
    
        G(old_parents,:) = 0;
        G(:,old_parents) = 0;        
        
    end
    
end