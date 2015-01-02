function get_points_omni_list(omni_file_dir, out_name, obj_list)

    if size(obj_list,1)>1
        obj_list = obj_list';
    end

    %get volume size
    fid = fopen([omni_file_dir '/projectMetadata.yaml'],'r');
%     fid = fopen(['./projectMetadata.yaml'],'r');
    
    metadata = fread(fid,'char')';
    
    data_dim_spot = strfind(metadata, 'dataDimensions');
    metadata = char(metadata(data_dim_spot(1):data_dim_spot(1)+100));
    % dataDimensions: [123, 3456, 12222]
    
    spaces = find(isspace(metadata));
    imSz(1) = str2double(metadata(spaces(1)+2:spaces(2)-2));
    imSz(2) = str2double(metadata(spaces(2)+1:spaces(3)-2));
    imSz(3) = str2double(metadata(spaces(3)+1:spaces(4)-2));
    
    
    chunk_size_spot = strfind(metadata, 'chunkDim');
    metadata = metadata(chunk_size_spot(1):end);
    spaces = find(isspace(metadata));
    chunkEdge = str2double(metadata(spaces(1)+1:spaces(2)-1));
    
    fclose(fid);
    
%
%

    MAX_CELL_IDS = length(obj_list);
    save_increment = 100000;
    save_buffer = 50000;

    
    max_chunk_size = chunkEdge*ones(1,3);
    
    num_chunks = ceil(imSz./max_chunk_size);
    c = zeros(1,3);
    
    
    fid = fopen([omni_file_dir '/segmentations/segmentation1/0/volume.uint32_t.raw']);
%     fid = fopen([omni_file_dir './volume.uint32_t.raw']);
    M = containers.Map(0,0);
    
    
    points = cell(MAX_CELL_IDS,1); %if we have more than 1000, will be slow
    num_points = zeros(MAX_CELL_IDS,1);
    num_saves = zeros(MAX_CELL_IDS,1);
    
    
    num_ids = 0;
    
    for n = 1:prod(num_chunks)
        tic
        
        [c(1) c(2) c(3)] = ind2sub(num_chunks, n);
        chunkSz = max_chunk_size;
        for k = 1:3
            if c(k) == num_chunks(k)
                chunkSz(k) = mod(imSz(k),max_chunk_size(k));
                if chunkSz(k) == 0;
                    chunkSz(k) = max_chunk_size(k);
                end
            end
        end
        
        chunk_data = fread(fid,prod(chunkSz),'uint32');
        
        chunk_data(chunk_data<=obj_min | chunk_data>=obj_max) = 0;
        
%         uni_ids = unique(chunk_data(:));
%         uni_ids(uni_ids==0) = [];
        
        
        
        for k = obj_list
            if ~M.isKey(k)
                num_ids = num_ids+1;
                M(k) = num_ids;
                points{num_ids} = zeros(save_increment+save_buffer,3);
                
            end
            
            id = M(k);
            
            p_ind = find(chunk_data==k);
            [px py pz] = ind2sub(chunkSz,p_ind);
            px = px + (c(1)-1)*max_chunk_size(1);
            py = py + (c(2)-1)*max_chunk_size(2);
            pz = pz + (c(3)-1)*max_chunk_size(3);
                                   
            points{id}(num_points(id)+1:num_points(id)+length(p_ind),:) = [px py pz];
            num_points(id) = num_points(id)+length(p_ind);
                        
            if num_points(id) >= save_increment                
                p = uint32(points{id}(1:num_points(id),:)); %#ok<NASGU>
                save(['./' out_name '_' num2str(k) '_' num2str(num_saves(id)) '.mat'],'p');
                num_points(id) = 0;
                points{id} = zeros(save_increment+save_buffer,3);
                num_saves(id) = num_saves(id)+1;

            end                    
            
            
        end
        
%         for k = 1:length(chunk_data)
%             if chunk_data(k)>=obj_min && chunk_data(k)<=obj_max
%                 d = chunk_data(k);
%                 if ~M.isKey(d)
%                     num_ids = num_ids+1;
%                     M(d) = num_ids;
%                     points{num_ids} = zeros(save_increment+save_buffer,3);
%                 endg
%                 id = M(d);
%                 
%                 pxyz = ind2sub(chunkSz,k) + (c-1).*max_chunk_size;
%                 points{id}(num_points(id),:) = pxyz;
%                 num_points(id) = num_points(id)+1;
%                 
%                 
%                 if num_points(id) >= save_increment                
%                     p = uint32(points{id}(1:num_points(id),:)); %#ok<NASGU>
%                     save(['./' out_name '_' num2str(d) '_' num2str(num_saves(id)) '.mat'],'p');
%                     num_points(id) = 0;
%                     points{id} = zeros(save_increment+save_buffer,3);
%                     num_saves(id) = num_saves(k)+1;
% 
%                 end       
%                 
%             end
%         end

        
        disp(['finishing block ' num2str(n) ' of ' num2str(prod(num_chunks))])
        toc
    end
    
    fclose(fid);
    
    
    ids = keys(M);
    id_vals = values(M);
    for c = 1:num_ids
        k = id_vals{c};            
        if k ~= 0
            p = uint32(points{k}(1:num_points(k),:));
            save(['./' out_name '_' num2str(ids{c}) '_' num2str(num_saves(k)) '.mat'],'p');                                  
        end
    end
    
    
end
    
    
    
    
    