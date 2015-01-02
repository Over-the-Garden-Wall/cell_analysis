function hull = make_locally_convex_hull(p, rad, remove_internal)

    old_hull = p(1,:);
    
    
    while 1
        
        hull = [];
        for n = 1:size(old_hull,1)
            
            d = sqrt((p(:,1)-old_hull(n,1)).^2 + (p(:,2)-old_hull(n,2)).^2);
            sub_p = p(d<=rad,:);
            
            if size(sub_p,1) > 2

                hull_inds  = convhull(sub_p(:,1), sub_p(:,2));
                sub_hull = sub_p(hull_inds,:);
                [sub_hull(:,1), sub_hull(:,2)] = poly2cw(sub_hull(:,1), sub_hull(:,2));

                if ~isempty(hull)
                    t_hull = [];
                    [t_hull(:,1), t_hull(:,2)] = polybool('union', hull(:,1), hull(:,2), sub_hull(:,1), sub_hull(:,2));
                    hull = t_hull;
                else
                    hull = sub_hull;        


                end

            end
        end
        
        if size(hull,1) == size(old_hull,1)
            if all(hull(:)==old_hull(:))
                break
            end
        end
        old_hull = hull;
    end
    
    if exist('remove_internal', 'var') && remove_internal
        nan_locs = find(isnan(hull(:,1)));

        while ~isempty(nan_locs)
            first_is_internal = inpolygon(hull(1,1), hull(1,2), hull(nan_locs(1)+1:end,1), hull(nan_locs(1)+1:end,2));
            if first_is_internal
                hull = hull(nan_locs(1)+1:end,:);
            elseif length(nan_locs)==1
                hull = hull(1:nan_locs(1)-1,:);
            else
                hull(nan_locs(1)+1:nan_locs(2),:) = [];
            end
            nan_locs = find(isnan(hull(:,1)));

        end
    end
    
            
end