function [Q P valid_coords phi psi] = find_planar_rotation_iterative(coords, allowed_outlier_percentage)    


    valid_coords = coords;
    
    valid_ts = false(size(coords,1),1);
    
    while sum(~valid_ts)>size(coords,1)*allowed_outlier_percentage

        P = fit_plane(valid_coords);
        P(1:3) = P(1:3)/sqrt(sum(P(1:3).^2));
        
        [Q phi psi] = find_planar_rotation(P);

        [dummy t] = project_onto_plane(valid_coords, P);

        valid_ts = abs(t-mean(t)) <= 2*std(t);

        valid_coords = valid_coords(valid_ts,:);
    
    end
    
    
end