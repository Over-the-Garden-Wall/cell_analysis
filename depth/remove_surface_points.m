function points = find_surface_points(cell_no)
    
    root_dir = '~/stratification/point_data/'; 

    fns = dir(root_dir);
    
    good_fns = [];
    num_points = 0;
    
    for n = 1:length(fns);
        fname = fns(n).name;
        if length(fname) > 6 && strcmp(fname(1:6),'points')
            under_loc = find(fname=='_');
            
            if str2double(fname(under_loc(end-1)+1:under_loc(end)-1))==cell_no
                good_fns(end+1) = n;
                load([root_dir fname]);
                num_points = num_points + size(p,1);
            end
        end
    end
    
    points = zeros(num_points,3);
    k = 0;
    for f = 1:length(good_fns)
        n = good_fns(f);
        load([root_dir fns(n).name]);
        points(k+1:k+size(p,1),:) = p;
        k = k+size(p,1);
    end
    
    % this shouldn't be necessary....
    points = unique(points,'rows');
    
            
    %sparse matrix plan...
    
%     fake_point_y = points(:,2)+points(:,3)*max_points(2);

    max_points = max(points)+2;

%     M = sparse(point_loc,ones(num_points,1), true(num_points,1), ...
%         max(points(:,3)+2)*2^32, 1);

    num_removed = 1;
    while num_removed > 0
    
        num_points = size(points,1);
        num_neighbors = zeros(num_points,1);

        M = sparse(points(:,1)+1,points(:,2)+1+points(:,3)*max_points(2), true(num_points,1), ...
            max_points(1)+2, max_points(2)*max_points(3));
        
        tic
        for n = 1:num_points % =\
            if rem(n,100000) == 0

                disp(['finished on #' num2str(n)]);
                toc
                tic
            end

           num_neighbors(n) = (...
                M(points(n,1), points(n,2)+points(n,3)*max_points(2)+1) + ...
                M(points(n,1)+2, points(n,2)+points(n,3)*max_points(2)+1) + ...
                M(points(n,1)+1, points(n,2)+points(n,3)*max_points(2)) + ...
                M(points(n,1)+1, points(n,2)+points(n,3)*max_points(2)+2) + ...
                M(points(n,1)+1, points(n,2)+(points(n,3)+1)*max_points(2)+1) + ...
                M(points(n,1)+1, points(n,2)+(points(n,3)-1)*max_points(2)+1));

        end

        to_remove = num_neighbors > 3 & num_neighbors < 6;
        num_removed = sum(to_remove);
        
        points = points(~to_remove,:);            
        disp(['removed ' num2str(num_removed) ' points.']);
        
    end
end
    