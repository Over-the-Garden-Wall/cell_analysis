function im = project_hdf5(fn, dim)
    
    a = hdf5info(fn);
    
    volSz = a.GroupHierarchy.Datasets.Dims;
    imSz = volSz;
    imSz(dim) = [];
    
    im = zeros(imSz);
    
    chunk_size = a.GroupHierarchy.Datasets.Chunksize;
    
    num_chunks = ceil(volSz./chunk_size);
    
    [chunk_nums(:,1), chunk_nums(:,2), chunk_nums(:,3)] = ind2sub(num_chunks,1:prod(num_chunks));
    
    for k = 1:size(chunk_nums,1)
        
        chunk = chunk_nums(k,:);
        chunk_bounds = [1, 1, 1; chunk_size] + ones(2,1)*((chunk-1).*chunk_size);
        chunk_bounds(2,:) = min([chunk_bounds(2,:); volSz]);
        
        minichunk = get_hdf5_file(fn, '/main', chunk_bounds(1,:), chunk_bounds(2,:));
        minichunk = squeeze(sum(minichunk, dim));

        im_bounds = chunk_bounds;
        im_bounds(:,dim) = [];
        im(im_bounds(1,1):im_bounds(2,1), im_bounds(1,2):im_bounds(2,2)) = im(im_bounds(1,1):im_bounds(2,1), im_bounds(1,2):im_bounds(2,2))+minichunk;
        
    end
end