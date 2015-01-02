function A = poly_area(poly_points)
    
    
    if any(isnan(poly_points))

%         warning('there are nans. Assuming they divide separate polygons. Internal polygons are treated as holes')
        
        nan_loc = find(isnan(poly_points(:,1)));
        nan_loc = [0; nan_loc; size(poly_points,1)+1];

        obj_size = nan_loc(2:end)-nan_loc(1:end-1);
        [dummy biggest_obj] = max(obj_size);
        
        num_polys = length(nan_loc)-1;
        p = cell(num_polys,1);
        first_p = zeros(num_polys,2);
        for n = 1:num_polys
            p{n} = poly_points(nan_loc(n)+1:nan_loc(n+1)-1,:);
            first_p(n,:) = p{n}(1,:);
        end
        
        is_internal = false(num_polys,1);
        
        for n = 1:num_polys
            is_internal = is_internal | (inpolygon(first_p(:,1), first_p(:,2), p{n}(:,1), p{n}(:,2)) & (1:num_polys)'~= n);
        end
        
        A = 0;
        for n = 1:num_polys
            if is_internal(n)
                A = A - poly_area(p{n});
            else
                A = A + poly_area(p{n});
            end            
        end
        
    
    else
        num_points = size(poly_points,1);
        for n = num_points:-1:2
            if all(poly_points(n,:)==poly_points(n-1,:))
                poly_points(n,:) = [];
            end
        end
        if ~(all(poly_points(end,:)==poly_points(1,:)))
            poly_points(end+1,:) = poly_points(1,:);
        end
        num_points = size(poly_points,1);
        
        
        if num_points <= 3
            A = 0;        
        else
        

            
            is_h = convhull(poly_points(:,1), poly_points(:,2));
            h = poly_points(is_h,:);

            A = poly_area_convex(h);        

            if size(h,1) ~= num_points
                is_in_h = false(num_points,1);
                is_in_h(is_h) = true;
                is_in_h(end) = is_in_h(1);

                num_comps = sum(is_in_h(1:end-1)~=is_in_h(2:end))/2;

                k_start = find(is_in_h,1,'first');
                new_order = [k_start:num_points-1 1:k_start];
                is_in_h = is_in_h(new_order);
                poly_points = poly_points(new_order,:);

                comp_begin = find(is_in_h(1:end-1) & ~is_in_h(2:end));
                comp_end = find(~is_in_h(1:end-1) & is_in_h(2:end))+1;

                for k = 1:num_comps
                    A = A - poly_area(poly_points(comp_begin(k):comp_end(k),:));
                end

            end

        end
    end
    
end



function A = poly_area_convex(poly_points)

    

    if ~all(poly_points(end,:)==poly_points(1,:))
        poly_points = [poly_points; poly_points(1,:)];
    end
    
    num_points = size(poly_points,1);
    
    A = 0;
    if num_points <= 3
        return
    end
    
%     m = mean(poly_points(1:end-1,:));
    m = poly_points(1,:);
    
    
    for n = 2:num_points-2
        a = sqrt(sum((m-poly_points(n,:)).^2));
        b = sqrt(sum((m-poly_points(n+1,:)).^2));
        c = sqrt(sum((poly_points(n+1,:)-poly_points(n,:)).^2));
        s = (a+b+c)/2;
        A = A + sqrt(s*(s-a)*(s-b)*(s-c));
    end
    
    if imag(A) ~= 0
%         warning('colinear points found in polygon...');
        A = real(A);
    end
    
end