function p = condense_points(p,rad)


    minp = min(p);
    for k = 1:size(p,2);
        p(:,k) = p(:,k) - minp(k);
        p(:,k) = 1+floor(p(:,k)/rad);
        
    end
    p = unique(p,'rows');
    
    for k = 1:size(p,2)
        p(:,k) = p(:,k)*rad + minp(k);
    end
end
    