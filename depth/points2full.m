function im = points2full(cell_no, dsmp_fact)
   
    root_dir = '~/stratification/point_data/'; 

    fns = dir(root_dir);
    
    good_fns = [];
%     num_points = 0;
    
    for n = 1:length(fns);
        fname = fns(n).name;
        if length(fname) > 6 && strcmp(fname(1:6),'points')
            under_loc = find(fname=='_');
            
            if str2double(fname(under_loc(end-1)+1:under_loc(end)-1))==cell_no
                good_fns(end+1) = n;
                load([root_dir fname]);
%                 num_points = num_points + size(p,1);
            end
        end
    end
    
    if ~isempty(good_fns)
        
        max_size = [0 0 0];
        for f = 1:length(good_fns)
            n = good_fns(f);
            load([root_dir fns(n).name]);
            max_size = max([max_size; max(p)]);
        end

        dsmp_fact = uint32(dsmp_fact);

        max_size = max_size./dsmp_fact+1;

        im = false(max_size);


        for f = 1:length(good_fns)
            n = good_fns(f);
            load([root_dir fns(n).name]);
            for k = 1:3
                p(:,k) = 1 + p(:,k)/dsmp_fact(k);
            end

            im(sub2ind(max_size,p(:,1),p(:,2),p(:,3))) = true;
        end
    
    else
        im = [];
    end
end
    
    

    