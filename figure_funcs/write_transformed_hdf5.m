function write_transformed_hdf5(cell_nos, out_resolution, out_fn);

    load('./T.mat');

    C = get_constants;
    
    root_dir = './point_data/';
    fns = dir(root_dir);
    is_valid_fn = false(length(fns),1);
    cell_nums = zeros(length(fns),1);
    
    %find valid fns
    for n = 1:length(fns)
        
        unders = strfind(fns(n).name,'_');
        if length(unders)<2
            cell_nums(n) = 0;
        else
            cell_nums(n) = str2double(fns(n).name(unders(1)+1:unders(2)-1));                
        end
        if any(cell_nums(n)==cell_nos)
            is_valid_fn(n) = true;
        end
        
    end
    
    valid_list = find(is_valid_fns');
    
   
    %find volume bounds
    
    min_val = Inf(1,3);
    max_val = -Inf(1,3);
    for k = cell_nos
        
        c_d = cell_data(k);
        p = c_d.get_surface;
        min_val = min([min_val; p]);
        max_val = max([max_val; p]);
        
    end
    
    vol_size = 1+ceil((max_val - min_val)./out_resolution);
    
    V = sparse([],[],[],prod(vol_size),1);
    
    for n = valid_list
        
        load([root_dir fns(n).name]);
        
        p = apply_transform(T,p);
        
        for k = 1:3
            p(:,k) = round((p(:,k)-min_val(k))/out_resolution(k))+1;
        end
        ind_list = sub2ind(vol_size, p(:,1), p(:,2), p(:,3));
        V(ind_list) = cell_nums(n);
        
    end
        