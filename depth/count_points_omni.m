function point_count = count_points_omni(omni_file_dir, obj_list)

    
%     fid = fopen([omni_file_dir '/segmentations/segmentation1/0/volume.uint32_t.raw']);
      page_file_dir = [omni_file_dir '/users/_default/segmentations/segmentation1/segments/'];
%       page_file_dir = './'

      fns = get_files_with_names_including(page_file_dir, 'data.ver4');
      
      max_ob = max(obj_list);
      num_obs = length(obj_list);
      point_count = zeros(num_obs,1,'uint32');
      
      id2ind = sparse(double(obj_list+1),ones(num_obs,1),(1:num_obs)',double(max_ob+1),1);
      
      for f = 1:length(fns)
         fid = fopen([page_file_dir '/' fns{f}]);
      
         f_data = fread(fid, 'uint32');
%          disp(fns{f})
         f_data = reshape(f_data, 12, length(f_data)/12);
         f_data = f_data';
             
         f_data(f_data(:,1) > max_ob,:) = [];
         inds = id2ind(f_data(:,1)+1);
         inds = full(inds);         
         f_data(inds==0,:) = [];
         inds(inds==0) = [];
         point_count(inds) = f_data(:,3);
%         for n = 1:size(f_data,1);
%             id = f_data(n,1);
%             if id < max_ob
%                 ind = id2ind(id+1);
%                 if ind > 0
%                     point_count(ind) = f_data(n,3);
%                 end
%             end
%         end
             
         
      
         fclose(fid);
            
            
      end

    
    
end
    
    
    
    
    