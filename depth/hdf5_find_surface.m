function hdf5_find_surface(fn, ids, out_name, chunk_size)

    save_increment = 100000;
    save_buffer = 50000;
    
    im_info = hdf5info(fn);
    imSz = im_info.GroupHierarchy.Datasets.Dims;
    
    
    num_chunks = ceil(imSz./chunk_size);
    
    
    points = cell(length(ids),1);
    num_points = zeros(length(ids),1);
    num_saves = zeros(length(ids),1);
    
    for n = 1:length(ids)
        points{n} = zeros(save_increment+save_buffer,3);
    end
    
    for n = 1:prod(num_chunks);
        
        [x y z] = ind2sub(num_chunks, n);
        
        start_coords = [x-1 y-1 z-1].*chunk_size;
        end_coords = [x y z].*chunk_size;
        end_coords = min([end_coords; imSz]);
        
        im = get_hdf5_file(fn,'/main', start_coords+1, end_coords);
    
        
        for k = 1:length(ids)

            p = find(im(:)==ids(k));
            [px py pz] = ind2sub(end_coords-start_coords,p);
            px = px + start_coords(1);
            py = py + start_coords(2);
            pz = pz + start_coords(3);

            points{k}(num_points(k)+1:num_points(k)+length(px),:) = [px py pz];
            num_points(k) = num_points(k)+length(px);
            


            if num_points(k) >= save_increment                
                p = uint32(points{k}(1:num_points(k),:));
                save(['./' out_name '_' num2str(k) '_' num2str(num_saves(k)) '.mat'],'p');
                num_points(k) = 0;
                points{k} = zeros(save_increment+save_buffer,3);
                num_saves(k) = num_saves(k)+1;

            end                                

        end
        
    disp(['finishing block ' num2str(n) ' of ' num2str(prod(num_chunks))])
    
    end

        
    
    for k = 1:length(ids)
        p = uint32(points{k}(1:num_points(k),:));
        save(['./' out_name '_' num2str(k) '_' num2str(num_saves(k)) '.mat'],'p');                                  
    end
    
    
end