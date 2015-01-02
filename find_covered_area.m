function [covered type_hulls] = find_covered_area(types)

    C = get_constants;
    
    type_hulls = cell(length(types),1);
    for k = 1:length(types);
        
        for c = C.type.(types{k})
            
            cell_dat = cell_data(c);
            
            hull = cell_dat.hull_2d;
            if isempty(type_hulls{k})
                type_hulls{k} = hull;
            else
                old_hull = type_hulls{k};
                type_hulls{k} = [];
                [type_hulls{k}(:,1) type_hulls{k}(:,2)] = polybool('union', hull(:,1), hull(:,2), old_hull(:,1), old_hull(:,2));
            end
            
        end
        
    end
    
    covered = type_hulls{1};
    for k = 2:length(types);
        was_covered = covered;
        covered = [];
        [covered(:,1) covered(:,2)] = polybool('intersection', was_covered(:,1), was_covered(:,2), type_hulls{k}(:,1), type_hulls{k}(:,2));
    end
        

end