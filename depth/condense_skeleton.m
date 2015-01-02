function [new_nodes edges] = condense_skeleton(nodes,edges, threshold)

    %dumb way, make better later; O(n^2)
    
    new_nodes = zeros(size(nodes));
    k = 0;
    
    max_node = size(nodes,1);
    is_valid = true(max_node,1);
    
    while any(is_valid)
        k = k+1;
        
        tn = nodes(find(is_valid,1,'first'),:);
        
        dist = sqrt((nodes(:,1)-tn(1)).^2 + ...
            (nodes(:,2)-tn(2)).^2 + ...
            (nodes(:,3)-tn(3)).^2);
        
        within_threshold = is_valid & dist <= threshold;
        new_nodes(k,1) = mean(nodes(within_threshold,1));
        new_nodes(k,2) = mean(nodes(within_threshold,2));
        new_nodes(k,3) = mean(nodes(within_threshold,3));
        
        for old_node = find(within_threshold')
            edges(edges==old_node) = max_node+k;
        end
        
        is_valid(within_threshold) = false;
%         nodes(within_threshold,:) = [];
    end
    
    edges = edges - max_node;        
    edges(edges(:,1)==edges(:,2),:) = [];
    
    
    new_nodes = new_nodes(1:k,:);
end