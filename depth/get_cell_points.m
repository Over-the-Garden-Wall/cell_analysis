function get_cell_points(cell_name, out_dir)


    if ~exist(['./' out_dir],'dir')
        mkdir(['./' out_dir]);
    end

    
    [seg_list_fns, omni_file_fns] = get_filenames_for_cell(cell_name);
    
    num_files = length(seg_list_fns);
    
    for n = 1:num_files
        
    disp(['working on ' seg_list_fns{n}]);
    tic
    
        valid_segs = read_concensus(seg_list_fns{n});
        get_points_omni_list(omni_file_fns{n}, ['./' out_dir '/' cell_name  '_' num2str(n) '_points'], valid_segs);
        [dends dendVals] = read_omni_dendrogram(omni_file_fns{n});
        save(['./' out_dir '/' cell_id '_' num2str(n) '_info.mat'],'valid_segs', 'dendVals', 'dends');

    toc
    
    end
    
end