function points = find_surface_points(cell_no, portion_untouchable)
    
    root_dir = '~/stratification/thinned_data/thin_'; 

    fn = [root_dir num2str(cell_no) '.mat'];
    
    load(fn)
            
    p = [1 0 -1; 0 1 -1; 1 -1 0];
    q = ones(3) - 2*eye(3);
    
    nhood = [eye(3); -eye(3); 1-eye(3); eye(3)-1; 1 1 1; -1 -1 -1; p; -p; q; -q];
    
    
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
                M(points(n,1)+1, points(n,2)+(points(n,3)-1)*max_points(2)+1)) + ;

        end

        to_remove = num_neighbors > 2 & num_neighbors < 6;
        num_removed = sum(to_remove);
        
        points = points(~to_remove,:);            
        disp(['removed ' num2str(num_removed) ' points.']);
        
    end
end
    