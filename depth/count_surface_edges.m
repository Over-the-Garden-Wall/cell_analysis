function surf_edge_count = count_surface_edges(cell_no)
    
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
    num_points = size(points,1);
    
            
    %sparse matrix plan...
    
%     fake_point_y = points(:,2)+points(:,3)*max_points(2);

    max_points = max(points)+2;

    M = sparse(points(:,1)+1,points(:,2)+1+points(:,3)*max_points(2), true(num_points,1), ...
        max_points(1)+2, max_points(2)*max_points(3));
%     M = sparse(point_loc,ones(num_points,1), true(num_points,1), ...
%         max(points(:,3)+2)*2^32, 1);

    surf_edge_count = 0;

    tic
    for n = 1:num_points % =\
        if rem(n,100000) == 0
            
            disp(['finished on #' num2str(n)]);
            toc
            tic
        end
        
        surf_edge_count = surf_edge_count + 6 - ...
            M(points(n,1), points(n,2)+points(n,3)*max_points(2)+1) - ...
            M(points(n,1)+2, points(n,2)+points(n,3)*max_points(2)+1) - ...
            M(points(n,1)+1, points(n,2)+points(n,3)*max_points(2)) - ...
            M(points(n,1)+1, points(n,2)+points(n,3)*max_points(2)+2) - ...
            M(points(n,1)+1, points(n,2)+(points(n,3)+1)*max_points(2)+1) - ...
            M(points(n,1)+1, points(n,2)+(points(n,3)-1)*max_points(2)+1);
            
    end
    
    
    
%    is_surface = ~(...
%         M(points(:,1), points(:,2)+points(:,3)*max_points(2)+1) & ...
%         M(points(:,1)+2, points(:,2)+points(:,3)*max_points(2)+1) & ...
%         M(points(:,1)+1, points(:,2)+points(:,3)*max_points(2)) & ...
%         M(points(:,1)+1, points(:,2)+points(:,3)*max_points(2)+2) & ...
%         M(points(:,1)+1, points(:,2)+(points(:,3)+1)*max_points(2)+1) & ...
%         M(points(:,1)+1, points(:,2)+(points(:,3)-1)*max_points(2)+1));

    
%    is_surface = ~(...
%         M(points(:,1), points(:,2)+points(:,3)*max_points(2)+1) & ...
%         M(points(:,1)+2, points(:,2)+points(:,3)*max_points(2)+1) & ...
%         M(points(:,1)+1, points(:,2)+points(:,3)*max_points(2)) & ...
%         M(points(:,1)+1, points(:,2)+points(:,3)*max_points(2)+2) & ...
%         M(points(:,1)+1, points(:,2)+(points(:,3)+1)*max_points(2)+1) & ...
%         M(points(:,1)+1, points(:,2)+(points(:,3)-1)*max_points(2)+1));

    
    
    
    % map plan failed, looks like checking for keys is really slow
%     point_loc = cell(num_points,1);
%     dummy = cell(num_points,1);
%     
%     for n = 1:num_points
%         %my god, matlab2009b doesn't support uint64 addition o.o
% %         point_loc{n} = uint64(coeff(1)*points(n,1)+1) + uint64(coeff(2)*points(n,2)) + uint64(coeff(3)*points(n,3));         
%         point_loc{n} = bitor(bitor(uint64(points(n,1)+1),bitshift(uint64(points(n,2)+1),16)), bitshift(uint64(points(n,3)+1),32));
%            
%         dummy{n} = true;
%     end
%     
%     
%     M = containers.Map(point_loc, dummy, 'uniformvalues', true); 
%     
% %     clear dummy
%     
%     is_surface = false(num_points,1);
%     
%     for n = 1:num_points
%         p = points(n,:);
%         
%         
%         
%         for k = 1:3
%             p(k) = p(k)+1;
%             p1 = bitor(bitor(uint64(p(1)+1),bitshift(uint64(p(2)+1),16)), bitshift(uint64(p(3)+1),32));
%             p(k) = p(k)-2;
%             p2 = bitor(bitor(uint64(p(1)+1),bitshift(uint64(p(2)+1),16)), bitshift(uint64(p(3)+1),32));
%             p(k) = p(k)+1;
%             
%             if ~M.isKey(p1) || ~M.isKey(p2)
%                 is_surface(n) = true;
%                 break
%             end
%         end
%     end 
    
end
    