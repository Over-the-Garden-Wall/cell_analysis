function write_transformed_hdf5(cell_nos, out_resolution, fn)

    try

    MAX_POINTS = 10^7;
    POINT_DUMP_THRESHOLD = MAX_POINTS/2;

    tic
    load('./T.mat');
    
    chunk_size = 128*ones(1,3);

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
    num_chunks = ceil(vol_size./chunk_size);
    
    all_p = zeros(MAX_POINTS, 4);    
    all_pk = 0;
    
    delete(fn);
    create_hdf5_file(fn, '/main', vol_size, chunk_size, [0 0 0]);
    
    for n = valid_list
        
%         n = valid_list(nk);
        
        load([root_dir fns(n).name]);
        
        p = apply_transform(T,double(p));
        
        for k = 1:3
            p(:,k) = round((p(:,k)-min_val(k))/out_resolution(k))+1;
        end
        
        p = unique(p,'rows');
        
        all_p(all_pk+ (1:size(p,1)),1:3) = p;
        all_p(all_pk+ (1:size(p,1)),4) = cell_nums(n);
        
        all_pk = all_pk+size(p,1);
        
        if all_pk > POINT_DUMP_THRESHOLD || n == valid_list(end)
           toc
           tic
            disp(['dumping ' num2str(all_pk) ' points']);
            
            all_p = all_p(1:all_pk,:);
            
            chunk_num = ceil(all_p(:,1)/chunk_size(1)) + ...
                num_chunks(1) * ceil(all_p(:,2)/chunk_size(2) - 1) + ...
                num_chunks(1)*num_chunks(2) * ceil(all_p(:,3)/chunk_size(3) - 1);
            
            t_list = unique(chunk_num);
%             t_list(t_list==0) = [];
            
            for t = t_list'
                sub_p = all_p(chunk_num==t,:);
                
                block_start = [0 0 0];
                [block_start(1), block_start(2), block_start(3)] = ind2sub(num_chunks, t);
                chunk_bounds = [1 1 1; chunk_size] + ones(2,1)*((block_start-1).*chunk_size);
                chunk_bounds(2,:) = min([chunk_bounds(2,:); vol_size]);
                
                for k = 1:3;
                    sub_p(:,k) = sub_p(:,k) - chunk_bounds(1,k) + 1;
                end
                
                chunk_vol = get_hdf5_file(fn, '/main', chunk_bounds(1,:), chunk_bounds(2,:));
                for k = 1:size(sub_p,1)
                    chunk_vol(sub_p(k,1), sub_p(k,2), sub_p(k,3)) = sub_p(k,4);
                end
                write_hdf5_file(fn, '/main', chunk_bounds(1,:), chunk_bounds(2,:), chunk_vol);
                
                
            end
            
            
            all_p = zeros(MAX_POINTS, 4);    
            all_pk = 0;
        end
        
    end
    
    catch ME
        disp(ME.message)
        disp('stop');
    end
    
end
        
    