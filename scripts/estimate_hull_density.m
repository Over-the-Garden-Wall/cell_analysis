figure; hold all
for k = 1:5; 
    h{k}=[]; 
    T{k} = 0; 
    cell_nums = C.type.(types{k}); 
    for c = cell_nums 
        cell_dat = cell_data(c); 
        hull = cell_dat.hull_2d;
        if isempty(h{k})
            h{k} = hull;
        else
            t = [];
            [t(:,1),t(:,2)] = polybool('union',h{k}(:,1), h{k}(:,2),hull(:,1),hull(:,2));
            h{k} = t;
        end
        T{k} = T{k} + poly_area(hull);
    end
    
%     h{k} = h{k}(~is_bad,:);
            
     plot(h{k}(:,1), h{k}(:,2));
    
    
    U{k} = poly_area(h{k});
    
end


