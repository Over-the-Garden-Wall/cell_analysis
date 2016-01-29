function write_transformed_hdf5(cell_nos, out_resolution, fn)

    try

    C = get_constants;
        
    MAX_POINTS = 10^7;
    POINT_DUMP_THRESHOLD = MAX_POINTS/2;

    tic
    load(C.trans_loc);
    
    chunk_size = 128*ones(1,3);

    
    root_dir = C.raw_point_dir;
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
    
    valid_list = find(is_valid_fn');
    
   
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
%     num_chunks = ceil(vol_size./chunk_size);
%     
%     all_p = zeros(MAX_POINTS, 4);    
%     all_pk = 0;
    
%     delete(fn);
%     create_hdf5_file(fn, '/main', vol_size, chunk_size, [0 0 0]);
    im = zeros(vol_size);

    for n = valid_list
        
%         n = valid_list(nk);
        
        load([root_dir fns(n).name]);
        
        p = apply_transform(T,double(p));
        
        for k = 1:3
            p(:,k) = round((p(:,k)-min_val(k))/out_resolution(k))+1;
        end
        
        p = unique(p,'rows');
        
        im(sub2ind(vol_size, p(:,1), p(:,2), p(:,3))) = cell_nums(n);        
        
    end
    hdf5write(fn, '/main', im);
    
    catch ME
        disp(ME.message)
        disp('stop');
    end
    
end
        
    